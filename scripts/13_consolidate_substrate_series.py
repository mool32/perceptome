#!/usr/bin/env python3
"""Phase 0 — consolidate the substrate-series (papers 4.3-4.8 + paper4 memory)
into `substrate_series_v1`: one long-format table of observed per-module effect
sizes, tagged by operation / context / effect-metric / verdict.

This is the calibration corpus for the Factor 2 (operation_intensity) estimator
of the perceptivity two-factor framework, and the regression-test corpus.

Deliberate scope limits (Phase 0):
  * Does NOT join HPA A_baseline yet — cell-type -> HPA mapping is a judgment
    call deferred to Phase 1 (column `A_baseline_hpa` left null as a placeholder).
  * Does NOT mix effect metrics — Cohen's d (4.3 immune detail, 4.5, paper4) and
    log2FC (4.4, 4.6-4.8) are tagged in `effect_metric`, never pooled.
  * 4.1 drugs / 4.2 cancer attractor are already consolidated elsewhere in the
    tool (pct.drugs, pct.reference.attractor) and use a different grain — they
    are referenced in STUDY_META but not ingested here.
"""
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/teo/Desktop/research/perceptual_modules")
TOOL = ROOT / "tool" / "perceptome2"
OUT = TOOL / "scripts"

# --- canonical module names from the v0.3 catalog (44 modules) ---
_cat = json.load(open(TOOL / "perceptome/catalog/data/modules_v03.json"))
CANON = set(_cat["modules"].keys())

# module-name reconciliation (paper4 memory uses older labels)
ALIAS = {"ERK/FOS": "ERK/MAPK", "NF-kB": "NF-κB", "NFkB": "NF-κB", "PI3K": "PI3K/PTEN"}


def canon(m):
    return ALIAS.get(str(m).strip(), str(m).strip())


def direction(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "na"
    if pd.isna(x):
        return "na"
    return "up" if x > 0.05 else "down" if x < -0.05 else "flat"


rows = []


def add(study, paradigm, operation, context, condition, cell_type, module,
        score_type, effect, metric, source, verdict=""):
    mc = canon(module)
    e = None
    if effect is not None and not (isinstance(effect, float) and pd.isna(effect)):
        e = round(float(effect), 4)
    rows.append(dict(
        study=study, paradigm=paradigm, operation=operation, context_class=context,
        condition=condition, cell_type=cell_type, module=mc,
        module_in_catalog=mc in CANON, score_type=score_type,
        effect=e, effect_metric=metric, direction=direction(effect),
        A_baseline_hpa=None,  # Phase-1 join placeholder
        verdict=verdict, source_file=source,
    ))


def rel(p):
    return str(Path(p).relative_to(ROOT))


# ---- paper4.3 immune activation: phase3 log2fc table (peak & 24h vs 0h) ----
f = ROOT / "paper4.3 immune/results/phase3/log2fc_table.csv"
df = pd.read_csv(f).rename(columns={"Unnamed: 0": "module"})
for _, r in df.iterrows():
    add("paper4.3", "immune_activation", "immune_activation", "active", "peak_vs_0h",
        "immune(activated)", r["module"], "activity", r["log2fc_peak_vs_0h"], "log2fc", rel(f))
    add("paper4.3", "immune_activation", "immune_activation", "active", "24h_vs_0h",
        "immune(activated)", r["module"], "activity", r["log2fc_24h_vs_0h"], "log2fc", rel(f))

# ---- paper4.4 muscle hypertrophy: RE activity log2fc per timepoint ----
f = ROOT / "paper4.4 muscle/results/02_re_activity_log2fc.csv"
for _, r in pd.read_csv(f).iterrows():
    add("paper4.4", "muscle_hypertrophy", "hypertrophy", "active", f"RE_{r['timepoint']}",
        "myofiber", r["module"], "activity", r["log2fc_mean"], "log2fc", rel(f))

# ---- paper4.6 organoid differentiation: active vs mature phase ----
f = ROOT / "paper4.6 organoid/results/02_phase_summary.csv"
df = pd.read_csv(f).rename(columns={"Unnamed: 0": "module"})
for _, r in df.iterrows():
    add("paper4.6", "organoid_differentiation", "differentiation", "active", "active_phase",
        "organoid", r["module"], "activity", r["active_lfc"], "log2fc", rel(f))
    add("paper4.6", "organoid_differentiation", "differentiation", "mature", "mature_phase",
        "organoid", r["module"], "activity", r["mature_lfc"], "log2fc", rel(f))

# ---- paper4.7 organoid regeneration: d1(active) vs d4(mature) ----
f = ROOT / "paper4.7 organoid_regen/results/02_full_module_log2fc.csv"
df = pd.read_csv(f).rename(columns={"Unnamed: 0": "module"})
for _, r in df.iterrows():
    add("paper4.7", "organoid_regeneration", "regeneration", "active", "d1_vs_d4",
        "organoid", r["module"], "activity", r["log2fc_d1_vs_d4"], "log2fc", rel(f))

# ---- paper4.8 organoid retinoid holdout: perturbations vs DMSO at d4 ----
f = ROOT / "paper4.8 organoid_holdout/results/02_log2fc_perturbation_vs_dmso.csv"
df = pd.read_csv(f).rename(columns={"Unnamed: 0": "module"})
for cond in ["atRA", "RXRi", "novita"]:
    for _, r in df.iterrows():
        add("paper4.8", "organoid_retinoid", "retinoid_perturbation", "drug_perturbation",
            cond, "organoid", r["module"], "activity", r[cond], "log2fc", rel(f))

# ---- paper4 memory consolidation: figure2 scorecard data (already long) ----
f = ROOT / "paper4 neuron perseption and memory/Memory Consolidation/repo/figures/figure2_scorecard_data.csv"
for _, r in pd.read_csv(f).iterrows():
    if str(r["module"]).upper() == "ALL":      # composite/orthogonality preds, not module effects
        continue
    add("paper4", "memory_consolidation", "memory_consolidation", "active",
        str(r["dataset"]), str(r["cell_type"]), r["module"], str(r["score_type"]),
        r["cohens_d"], "cohens_d", rel(f), verdict=str(r["verdict"]))

# ---- paper4.5 intestinal homeostasis: single-module d-values from scorecard ----
f = ROOT / "paper4.5 epithelial/results/02_scorecard.csv"
pat = re.compile(r"^(?P<cell>\w[\w ]*?) vs ISC: (?P<mod>[\w\-/κβα]+) d")
for _, r in pd.read_csv(f).iterrows():
    m = pat.search(str(r["observation"]))
    if m and pd.notna(r["value"]):
        add("paper4.5", "intestinal_homeostasis", "terminal_differentiation", "terminal/homeostasis",
            f"{m.group('cell')}_vs_ISC", m.group("cell"), m.group("mod"), "activity",
            r["value"], "cohens_d", rel(f), verdict=str(r["verdict"]))

obs = pd.DataFrame(rows)

# ---- study-level metadata + capacity-floor verdict (primary pass rate from scorecards) ----
STUDY_META = {
    "paper4":   ("memory consolidation (hippocampus)", "active differentiation", "PASS (activity layer)"),
    "paper4.3": ("immune activation (ILC3/general)",   "active stimulation",      "FAIL (sentinel not reproduced)"),
    "paper4.4": ("muscle resistance exercise",         "active hypertrophy",      "PASS 5/6 (peak-aligned)"),
    "paper4.5": ("intestinal epithelium (adult)",      "terminal/homeostasis",    "FAIL 1/8 (BOUNDARY)"),
    "paper4.6": ("intestinal organoid differentiation","active differentiation",  "FAIL 2/7 (null validated)"),
    "paper4.7": ("intestinal organoid regeneration",   "active regeneration",     "STRONG PASS 5/5"),
    "paper4.8": ("organoid retinoid holdout",          "drug perturbation",       "PARTIAL 3/5"),
}

# ---- write artifacts ----
obs.to_csv(OUT / "substrate_series_v1.csv", index=False)
payload = {
    "schema_version": "substrate_series_v1",
    "grain": "one row per (study, condition, cell_type, module, score_type) observed effect",
    "metrics": "effect_metric tags d vs log2fc; NEVER pooled across metrics",
    "deferred": "A_baseline_hpa join (cell->HPA mapping) deferred to Phase 1",
    "study_meta": {k: dict(zip(["paradigm", "context_class", "capacity_floor_verdict"], v))
                   for k, v in STUDY_META.items()},
    "n_observations": int(len(obs)),
    "observations": obs.where(pd.notna(obs), None).to_dict(orient="records"),
}
json.dump(payload, open(OUT / "substrate_series_v1.json", "w"), ensure_ascii=False, indent=1)

# ================= REVIEW PRINTOUT =================
pd.set_option("display.width", 240)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_colwidth", 40)

print(f"\nTotal observations: {len(obs)}  |  written to substrate_series_v1.{{csv,json}}")

print("\n=== STUDY SUMMARY ===")
summ = (obs.groupby(["study", "operation", "context_class", "effect_metric"])
        .agg(n_obs=("effect", "size"), n_modules=("module", "nunique"),
             n_nonnull=("effect", "count")).reset_index())
summ["paradigm"] = summ["study"].map(lambda s: STUDY_META[s][0])
summ["capacity_verdict"] = summ["study"].map(lambda s: STUDY_META[s][2])
print(summ.to_string(index=False))

print("\n=== MODULE-NAME COVERAGE CHECK ===")
bad = sorted(obs.loc[~obs["module_in_catalog"], "module"].unique())
print("modules NOT in v0.3 catalog:", bad if bad else "none — all reconciled")

print("\n=== FACTOR-2 ANCHOR VIEW: infrastructure/capacity modules across operations ===")
focus = ["HSF1", "UPR-ATF6", "UPR-IRE1", "UPR-PERK", "NRF2", "NFAT", "mTOR", "ERK/MAPK", "Cell Cycle"]
av = obs[obs["module"].isin(focus)].copy()
# one representative active-context value per (study, condition, module)
piv = (av.pivot_table(index=["study", "operation", "condition", "effect_metric"],
                      columns="module", values="effect", aggfunc="mean")
       .reindex(columns=focus))
print(piv.round(3).to_string())
