"""EXPLORATORY (not pre-registered, not locked): dose-response and buffering in Replogle K562_gwps.

These data were already used by the locked loop-integrity test, so nothing here
is confirmatory. The purpose is hypothesis generation: anything that looks
interesting must be re-tested under a new pre-registration on an independent screen.

Part A  readability   — can the residual level of each knocked-down gene be read per cell?
Part B  dose-response — within a knockdown population, does the module output follow the
                        per-cell residual level of the knocked-down gene?
Part C  buffering     — across all ~9.8k knockdowns (pseudobulk): how much does the
                        transcriptome move per knockdown, and how local is the response to
                        the module the gene belongs to? Summarised by perceptome module and
                        by the Paper 7 evolutionary cohorts.

Usage:
  python scripts/explore/15_dose_response_buffering.py \
      --sc   ~/data/replogle2022/K562_gwps_raw_singlecell_01.h5ad \
      --bulk ~/data/replogle2022/K562_gwps_normalized_bulk_01.h5ad \
      --out  results/explore_dose_buffering
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import perceptome as pct  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CONTRASTS = ROOT / "data" / "evidence" / "replogle2022_K562_gwps_contrasts.tsv"
SEED = 0
N_CTRL = 5000
CONTROL = "non-targeting"

# Dose-response cases: (knockdown, module whose output is read). The first block are the
# clearest known breaks/inputs from the locked run; the rest are every loop knockdown.
CASES_FIRST = [("KEAP1", "NRF2"), ("CUL3", "NRF2"), ("NFE2L2", "NRF2"), ("HSPA5", "UPR-ATF6"),
               ("ATF6", "UPR-ATF6"), ("SCAP", "SREBP"), ("INSIG1", "SREBP"), ("MBTPS2", "SREBP"),
               ("VHL", "HIF"), ("HIF1A", "HIF"), ("CISH", "JAK-STAT"), ("STAT5B", "JAK-STAT")]

# Paper 7 (evolution) doubly-robust cohorts, as reported in SESSION_2026-05-15_summary.
EVO_COHORTS = {
    "ancient_unit": ["HSF1", "ERK/MAPK", "PI3K/PTEN", "AMPK"],
    "young_anti_unit": ["p53", "VDR", "cGAS-STING", "NFAT"],
}


def machinery(module):
    """Genes whose knockdown counts as perturbing the module (not the huge cascade lists)."""
    g = set()
    for s in ("core", "sensor", "tf", "feedback"):
        g |= set(pct.get_genes(module, s))
    g |= set(pct.get_loop(module)["feedback_core"])
    return sorted(g)


def write_tsv(path, rows, cols):
    with open(path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


# ---------------- single cell: parts A and B ----------------
def load_sc(path):
    import anndata as ad
    import h5py
    a = ad.read_h5ad(path, backed="r")
    labels = a.obs["gene"].astype(str).to_numpy()
    genes = a.var["gene_name"].astype(str).to_numpy()
    h = h5py.File(path, "r")
    return labels, genes, h["X"]


def read_rows(X, idx, chunk=2000):
    idx = np.sort(np.asarray(idx))
    out = [np.asarray(X[idx[i:i + chunk]]) for i in range(0, len(idx), chunk)]
    return idx, np.vstack(out).astype(np.float32)


def lognorm(C):
    tot = C.sum(1, keepdims=True)
    return np.log1p(C / np.maximum(tot, 1) * 1e4), np.log(np.maximum(tot[:, 0], 1))


def resid(y, x):
    b = np.polyfit(x, y, 1)
    return y - np.polyval(b, x)


def parts_ab(sc_path, out):
    labels, genes, X = load_sc(sc_path)
    gpos = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(SEED)
    ctrl = np.sort(rng.choice(np.where(labels == CONTROL)[0], N_CTRL, replace=False))
    _, Cc = read_rows(X, ctrl)
    Lc, ltc = lognorm(Cc)
    mu, sd = Lc.mean(0), Lc.std(0)
    sd[sd == 0] = 1

    loop_rows = list(csv.DictReader(open(CONTRASTS), delimiter="\t"))
    cases = list(CASES_FIRST)
    for r in loop_rows:
        if r["in_screen"] == "yes" and (r["gene"], r["module"]) not in cases:
            cases.append((r["gene"], r["module"]))

    A, B, BINS = [], [], []
    for kd, mod in cases:
        cells = np.where(labels == kd)[0]
        if kd not in gpos or len(cells) < 30:
            A.append({"kd": kd, "module": mod, "n_cells": len(cells),
                      "kd_gene_measured": kd in gpos, "readable": False})
            continue
        _, Ck = read_rows(X, cells)
        Lk, ltk = lognorm(Ck)
        g = gpos[kd]
        ctrl_raw = Cc[:, g]
        a_row = {"kd": kd, "module": mod, "n_cells": len(cells), "kd_gene_measured": True,
                 "ctrl_mean_umi": float(ctrl_raw.mean()), "ctrl_frac_nonzero": float((ctrl_raw > 0).mean()),
                 "kd_frac_nonzero": float((Ck[:, g] > 0).mean())}
        a_row["readable"] = bool(a_row["ctrl_mean_umi"] >= 1.0 and a_row["ctrl_frac_nonzero"] >= 0.5)
        A.append(a_row)
        outg = [gpos[x] for x in pct.loop_output_genes(mod) if x in gpos and x != kd]
        if len(outg) < 2:
            continue
        score_k = ((Lk[:, outg] - mu[outg]) / sd[outg]).mean(1)
        score_c = ((Lc[:, outg] - mu[outg]) / sd[outg]).mean(1)
        xk, xc = Lk[:, g], Lc[:, g]
        row = {"kd": kd, "module": mod, "readable": a_row["readable"], "n_cells": len(cells),
               "mean_score_kd": float(score_k.mean()),
               "rho_kd": float(stats.spearmanr(xk, score_k).statistic),
               "rho_ctrl": float(stats.spearmanr(xc, score_c).statistic),
               "rho_kd_partial_depth": float(stats.spearmanr(resid(xk, ltk), resid(score_k, ltk)).statistic),
               "rho_ctrl_partial_depth": float(stats.spearmanr(resid(xc, ltc), resid(score_c, ltc)).statistic)}
        B.append(row)
        # bins: zero, then tertiles of nonzero residual level; controls as reference
        nz = xk > 0
        edges = np.quantile(xk[nz], [0, 1 / 3, 2 / 3, 1]) if nz.sum() >= 15 else None
        BINS.append({"kd": kd, "module": mod, "bin": "controls", "n": len(score_c),
                     "mean_kd_gene": float(xc.mean()), "mean_score": float(score_c.mean()),
                     "se_score": float(score_c.std() / np.sqrt(len(score_c)))})
        groups = [("zero", ~nz)]
        if edges is not None:
            for i, nm in enumerate(["low", "mid", "high"]):
                lo, hi = edges[i], edges[i + 1]
                groups.append((nm, nz & (xk >= lo) & ((xk <= hi) if i == 2 else (xk < hi))))
        for nm, m in groups:
            if m.sum() == 0:
                continue
            BINS.append({"kd": kd, "module": mod, "bin": nm, "n": int(m.sum()),
                         "mean_kd_gene": float(xk[m].mean()), "mean_score": float(score_k[m].mean()),
                         "se_score": float(score_k[m].std() / np.sqrt(m.sum()))})

    write_tsv(out / "A_readability.tsv", A, ["kd", "module", "n_cells", "kd_gene_measured", "ctrl_mean_umi",
                                             "ctrl_frac_nonzero", "kd_frac_nonzero", "readable"])
    write_tsv(out / "B_dose_response.tsv", B, ["kd", "module", "readable", "n_cells", "mean_score_kd", "rho_kd",
                                               "rho_ctrl", "rho_kd_partial_depth", "rho_ctrl_partial_depth"])
    write_tsv(out / "B_dose_bins.tsv", BINS, ["kd", "module", "bin", "n", "mean_kd_gene", "mean_score", "se_score"])
    return A, B


# ---------------- pseudobulk: part C ----------------
def part_c(bulk_path, out):
    import anndata as ad
    a = ad.read_h5ad(bulk_path)          # pseudobulk is small; X = z-normalised expression per perturbation
    genes = (a.var["gene_name"].astype(str).to_numpy() if "gene_name" in a.var
             else np.asarray(a.var_names, dtype=str))
    gpos = {g: i for i, g in enumerate(genes)}
    lab = np.asarray(a.obs_names, dtype=str)
    sym = np.array([s.split("_")[1] if "_" in s else s for s in lab])
    ncell = a.obs["num_cells_filtered"].to_numpy()
    fold = a.obs["fold_expr"].to_numpy() if "fold_expr" in a.obs else np.full(len(lab), np.nan)
    Z = np.asarray(a.X, dtype=np.float32)
    base = a.var["mean"].to_numpy() if "mean" in a.var else np.full(len(genes), np.nan)

    # one row per gene: the guide group with most cells
    best = {}
    for i, s in enumerate(sym):
        if s == CONTROL:
            continue
        if s not in best or ncell[i] > ncell[best[s]]:
            best[s] = i

    mods = pct.list_modules()
    act = {m: [gpos[g] for g in pct.get_genes(m, "activity") if g in gpos] for m in mods}
    owner = {}
    for m in mods:
        for g in machinery(m):
            owner.setdefault(g, []).append(m)

    rows = []
    for s, i in best.items():
        z = Z[i].copy()
        if s in gpos:
            z[gpos[s]] = 0.0
        rms = float(np.sqrt(np.mean(z ** 2)))
        r = {"kd": s, "n_cells": int(ncell[i]), "fold_expr": float(fold[i]), "rms_response": rms,
             "modules": ";".join(owner.get(s, []))}
        for m in owner.get(s, []):
            idx = [j for j in act[m] if genes[j] != s]
            if idx:
                own = float(np.mean(z[idx] ** 2))
                r[f"locality::{m}"] = own / max(np.mean(z ** 2), 1e-12)
        rows.append(r)

    eff = [r for r in rows if r["n_cells"] >= 50 and np.isfinite(r["fold_expr"]) and r["fold_expr"] <= 0.5]
    allrms = np.array([r["rms_response"] for r in eff])
    for r in eff:
        r["rms_pctile_among_efficient"] = float((allrms < r["rms_response"]).mean())
    per_mod = []
    for m in mods:
        kds = [r for r in eff if m in r["modules"].split(";")]
        mb = [base[j] for j in act[m]]
        per_mod.append({
            "module": m, "n_kd_efficient": len(kds),
            "kds": ",".join(r["kd"] for r in kds),
            "median_rms_pctile": float(np.median([r["rms_pctile_among_efficient"] for r in kds])) if kds else "",
            "median_locality": float(np.median([r[f"locality::{m}"] for r in kds if f"locality::{m}" in r]))
            if any(f"locality::{m}" in r for r in kds) else "",
            "baseline_activity_gene_mean": float(np.nanmean(mb)) if mb else "",
            "evo_cohort": next((k for k, v in EVO_COHORTS.items() if m in v), ""),
        })
    write_tsv(out / "C_buffering_per_kd.tsv", eff, ["kd", "n_cells", "fold_expr", "rms_response",
                                                     "rms_pctile_among_efficient", "modules"])
    write_tsv(out / "C_buffering_per_module.tsv", per_mod, ["module", "evo_cohort", "n_kd_efficient", "kds",
                                                           "median_rms_pctile", "median_locality",
                                                           "baseline_activity_gene_mean"])
    return per_mod, len(eff)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sc", required=True)
    p.add_argument("--bulk", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    per_mod, n_eff = part_c(a.bulk, out)
    A, B = parts_ab(a.sc, out)
    summary = {
        "status": "EXPLORATORY — hypothesis generation only, data already used by a locked test",
        "A_n_cases": len(A), "A_n_readable": sum(1 for r in A if r.get("readable")),
        "B_n_cases": len(B), "C_n_efficient_kds": n_eff,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
