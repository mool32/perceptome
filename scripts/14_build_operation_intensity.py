#!/usr/bin/env python3
"""Phase 1 — build the Factor 2 reference `operation_intensity_v1` from
substrate_series_v1.

operation_intensity(op, module) := the PEAK engagement (ramp) that operation
`op` drives in `module`, measured in calibration cells that have capacity.
This is the multiplicand against a new cell's Factor-1 capacity:

    predicted_engagement(new_cell, op, module)
        = capacity_floor(new_cell, module)  ×  operation_intensity(op, module)

Honest scope (locked): data-limited (n≈10 substrate paradigms, strict 3a FAILed).
Shipped as ORDINAL intensity bands + direction + n_obs, NEVER a point estimate.
Metrics (log2fc vs cohens_d) are kept separate — never pooled.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/teo/Desktop/research/perceptual_modules")
TOOL = ROOT / "tool" / "perceptome2"
OUT = TOOL / "scripts"

import sys
sys.path.insert(0, str(TOOL))
import perceptome as pct  # noqa: E402

# ---- substrate cell -> HPA cell type mapping (for A_baseline QC of calibration cells) ----
# intestinal cells map exactly; muscle/immune/neuron map by lineage; organoid is a
# heterogeneous proliferative mix -> enteric stem compartment as proxy (flagged).
CELL_TO_HPA = {
    "myofiber": "myonuclei",
    "organoid": "enteric stem cells",
    "immune(activated)": "innate lymphoid cells",
    "DG": "brain excitatory neurons",
    "CA1": "brain excitatory neurons",
    "all_hippocampus": "brain excitatory neurons",
    "Goblet": "goblet cells",
    "Enterocyte": "enterocytes",
    "Paneth": "paneth cells",
}
MAP_CONFIDENCE = {
    "myofiber": "high", "organoid": "medium (heterogeneous→proliferative proxy)",
    "immune(activated)": "high (ILC3→ILC)", "DG": "high", "CA1": "high",
    "all_hippocampus": "medium", "Goblet": "exact", "Enterocyte": "exact", "Paneth": "exact",
}

SAT_HI = 4.5  # capacity-floor saturation threshold (floor.py default)


def band(absval):
    """Ordinal intensity band from |effect| (within-metric)."""
    if absval is None or np.isnan(absval):
        return "no_data"
    if absval < 0.10:
        return "none"
    if absval < 0.30:
        return "low"
    if absval < 0.60:
        return "medium"
    if absval < 1.00:
        return "high"
    return "extreme"


obs = pd.DataFrame(json.load(open(OUT / "substrate_series_v1.json"))["observations"])
obs = obs[obs["module_in_catalog"] & obs["effect"].notna()].copy()  # catalog modules only
obs["abs"] = obs["effect"].abs()

# ---- A_baseline of the calibration cell (QC: was the cell already saturated?) ----
ref = pct.load_hpa_perceptivity()
A = ref["A"]


def calib_Abaseline(row):
    hpa = CELL_TO_HPA.get(row["cell_type"])
    if hpa is None or hpa not in A.index or row["module"] not in A.columns:
        return np.nan
    return float(A.loc[hpa, row["module"]])


obs["calib_A_baseline"] = obs.apply(calib_Abaseline, axis=1)

# ---- collapse to peak engagement per (operation, module, metric) ----
recs = []
for (op, mod, metric), g in obs.groupby(["operation", "module", "effect_metric"]):
    i = g["abs"].idxmax()                       # peak |effect| row
    peak = float(g.loc[i, "effect"])
    cond = g.loc[i, "condition"]
    cell = g.loc[i, "cell_type"]
    calA = g.loc[i, "calib_A_baseline"]
    direction = "suppressive" if peak <= -0.10 else ("up" if peak >= 0.10 else "flat")
    recs.append(dict(
        operation=op, module=mod, effect_metric=metric,
        peak_effect=round(peak, 4), intensity_band=band(abs(peak)), direction=direction,
        n_obs=int(len(g)), peak_condition=cond, calib_cell=cell,
        calib_A_baseline=None if np.isnan(calA) else round(calA, 3),
        calib_saturated=bool(calA > SAT_HI) if not np.isnan(calA) else None,
    ))
ref_df = pd.DataFrame(recs).sort_values(["operation", "intensity_band", "module"])

# ---- write Factor-2 reference ----
payload = {
    "schema_version": "operation_intensity_v1",
    "definition": "peak engagement (ramp) operation drives in module, in capacity-bearing calibration cells",
    "bands": {"none": "<0.10", "low": "0.10-0.30", "medium": "0.30-0.60", "high": "0.60-1.0", "extreme": ">1.0"},
    "metrics_separate": True,
    "disclaimer": "calibration-grade, not statistically validated (n≈10 paradigms, strict 3a FAILed). Ordinal only.",
    "cell_to_hpa_map": CELL_TO_HPA,
    "operations": sorted(obs["operation"].unique().tolist()),
    "table": ref_df.where(pd.notna(ref_df), None).to_dict(orient="records"),
}
json.dump(payload, open(OUT / "operation_intensity_v1.json", "w"), ensure_ascii=False, indent=1)
ref_df.to_csv(OUT / "operation_intensity_v1.csv", index=False)
# canonical bundled copy read by perceptome.load_operation_intensity()
PKG_DATA = TOOL / "perceptome" / "perceptivity" / "data"
json.dump(payload, open(PKG_DATA / "operation_intensity_v1.json", "w"), ensure_ascii=False, indent=1)

# ================= REVIEW PRINTOUT =================
pd.set_option("display.width", 240); pd.set_option("display.max_columns", 20); pd.set_option("display.max_colwidth", 34)

print("\n=== A_baseline → HPA mapping (calibration cells) ===")
mp = pd.DataFrame([dict(substrate_cell=k, hpa_cell=v, confidence=MAP_CONFIDENCE[k]) for k, v in CELL_TO_HPA.items()])
print(mp.to_string(index=False))

print(f"\n=== operation_intensity_v1: {len(ref_df)} (operation × module × metric) entries ===")
print("band counts:", ref_df["intensity_band"].value_counts().to_dict())
print("direction counts:", ref_df["direction"].value_counts().to_dict())

print("\n=== operation-level intensity profile (focus modules; peak |effect|) ===")
focus = ["HSF1", "UPR-ATF6", "UPR-PERK", "UPR-IRE1", "NRF2", "NFAT", "mTOR", "ERK/MAPK", "Cell Cycle", "Wnt"]
piv = (ref_df[ref_df["module"].isin(focus)]
       .pivot_table(index=["operation", "effect_metric"], columns="module", values="peak_effect", aggfunc="first")
       .reindex(columns=focus))
print(piv.round(2).to_string())

print("\n=== suppressive (down) intensities — the upward-asymmetric cases ===")
sup = ref_df[ref_df["direction"] == "suppressive"].sort_values("peak_effect")
print(sup[["operation", "module", "peak_effect", "intensity_band", "peak_condition"]].head(15).to_string(index=False))

print("\n=== QC: calibration cells flagged saturated (intensity may be ceiling-limited) ===")
sat = ref_df[ref_df["calib_saturated"] == True]
print(f"{len(sat)} of {len(ref_df)} entries from saturated calibration cells")
if len(sat):
    print(sat[["operation", "module", "calib_cell", "calib_A_baseline", "peak_effect"]].head(12).to_string(index=False))
