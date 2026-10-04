"""Candidate inventory for a confirmatory dose-response test in Replogle RPE1 (metadata only).

Reads:
  - single-cell obs/var (knockdown labels, gem groups, cells per knockdown, readout genes);
  - pseudobulk obs (author knockdown efficiency fold_expr);
  - raw counts of the knocked-down genes ONLY, in up to 5,000 non-targeting control cells,
    to judge per-cell readability. No module output gene, no knockdown cell and no
    knockdown-vs-output relation is read: the hypothesis to be tested stays unseen.

Usage:
  python scripts/explore/16_rpe1_candidates.py \
      --sc   ~/data/replogle2022/rpe1_raw_singlecell_01.h5ad \
      --bulk ~/data/replogle2022/rpe1_normalized_bulk_01.h5ad \
      --out  results/rpe1_candidates
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import perceptome as pct  # noqa: E402

SEED = 0
N_CTRL = 5000
CONTROL = "non-targeting"
MIN_CELLS = 50
MAX_FOLD = 0.50
READ_UMI = 1.0           # readable: control mean >= 1 UMI per cell ...
READ_FRAC = 0.5          # ... and nonzero in >= 50 % of control cells (same rule as the K562 exploration)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sc", required=True)
    p.add_argument("--bulk", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    import h5py
    sc = ad.read_h5ad(a.sc, backed="r")
    obs = sc.obs
    gcol = "gene" if "gene" in obs else "perturbation"
    labels = obs[gcol].astype(str).to_numpy()
    genes = (sc.var["gene_name"].astype(str).to_numpy() if "gene_name" in sc.var
             else np.asarray(sc.var_names, dtype=str))
    gpos = {g: i for i, g in enumerate(genes)}
    ncell = dict(zip(*np.unique(labels, return_counts=True)))

    bulk = ad.read_h5ad(a.bulk, backed="r")
    bsym = np.array([s.split("_")[1] if "_" in s else s for s in np.asarray(bulk.obs_names, dtype=str)])
    bn = bulk.obs["num_cells_filtered"].to_numpy()
    bf = bulk.obs["fold_expr"].to_numpy()
    fold = {}
    for s, n, f in zip(bsym, bn, bf):          # guide group with most cells
        if s not in fold or n > fold[s][0]:
            fold[s] = (n, f)

    # every perturbation gene of every loop, with its role
    cand = []
    for m in pct.list_modules():
        e = pct.get_loop(m)
        for cls, gl in e["perturbations"].items():
            for g in gl:
                cand.append((g, m, cls))
    kd_genes = sorted({g for g, _, _ in cand if g in gpos})

    # control-cell counts of the knocked-down genes only
    rng = np.random.default_rng(SEED)
    ctrl = np.where(labels == CONTROL)[0]
    ctrl = np.sort(rng.choice(ctrl, min(N_CTRL, len(ctrl)), replace=False))
    cols = np.array(sorted(gpos[g] for g in kd_genes))
    h = h5py.File(a.sc, "r")
    X = h["X"]
    if not isinstance(X, h5py.Dataset):
        raise SystemExit("X is not a dense dataset; stop and report (the K562 file was dense)")
    sub = np.vstack([np.asarray(X[ctrl[i:i + 1000]])[:, cols] for i in range(0, len(ctrl), 1000)])
    col_of = {int(c): k for k, c in enumerate(cols)}

    rows = []
    for g, m, cls in cand:
        r = {"module": m, "class": cls, "gene": g, "n_cells": int(ncell.get(g, 0)),
             "fold_expr": float(fold[g][1]) if g in fold else "", "kd_gene_measured": g in gpos}
        if g in gpos:
            v = sub[:, col_of[gpos[g]]]
            r["ctrl_mean_umi"] = float(v.mean())
            r["ctrl_frac_nonzero"] = float((v > 0).mean())
            r["readable"] = bool(r["ctrl_mean_umi"] >= READ_UMI and r["ctrl_frac_nonzero"] >= READ_FRAC)
        else:
            r["readable"] = False
        outg = [x for x in pct.loop_output_genes(m) if x in gpos and x != g]
        r["n_output_measured"] = len(outg)
        r["eligible"] = bool(r["n_cells"] >= MIN_CELLS and r["fold_expr"] != "" and r["fold_expr"] <= MAX_FOLD
                             and r["readable"] and len(outg) >= 2)
        rows.append(r)

    cols_out = ["module", "class", "gene", "n_cells", "fold_expr", "kd_gene_measured", "ctrl_mean_umi",
                "ctrl_frac_nonzero", "readable", "n_output_measured", "eligible"]
    with open(out / "rpe1_candidates.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols_out)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols_out])
    summary = {
        "obs_columns": list(map(str, obs.columns)),
        "knockdown_column": gcol,
        "n_cells_total": int(len(labels)),
        "n_control_cells": int((labels == CONTROL).sum()),
        "n_readout_genes": int(len(genes)),
        "n_candidate_rows": len(rows),
        "n_present_in_screen": sum(1 for r in rows if r["n_cells"] > 0),
        "n_readable": sum(1 for r in rows if r["readable"]),
        "eligible": sorted({f"{r['gene']}|{r['module']}|{r['class']}" for r in rows if r["eligible"]}),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
