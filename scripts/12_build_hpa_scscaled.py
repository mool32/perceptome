"""Build single-cell-scaled HPA perceptivity reference (perceptome v0.2.3).

The bundled hpa_perceptivity_v03.npz reference is in **pseudobulk units** —
each (cell type, module) value is mean log1p(nCPM) over the module's core/
activity genes computed directly from HPA's per-cell-type pseudobulk TSV.

This is the right scale for cell-type-level analyses. But when the user wants
to project a per-cell single-cell tumor dataset into the eigenspace
(score_modules + z-score against HPA + project), the pseudobulk-scale
reference is mismatched: per-cell scores after the standard
normalize_total(1e4) + log1p preprocessing fall in a 0–0.5 range, while the
pseudobulk HPA reference sits in a 1–3 range. Z-scoring per-cell tumor data
against pseudobulk μ, σ produces extremely negative z-vectors that collapse
to the "low-engagement" pole of the eigenspace (top-1 nearest HPA cell type
becomes erythrocytes for all cells — a calibration artifact).

This script generates a complementary reference where HPA pseudobulk nCPM is
**passed through the same normalize_total(1e4) + log1p preprocessing as
single-cell tumor data**, then module-scored. The resulting μ, σ are on the
same scale as per-cell tumor scores, and per-cell projection works correctly.

Outputs (under perceptome/perceptivity/data/):
  hpa_perceptivity_v03_scscaled.npz   R, A, C, headroom (154 × 44)
  hpa_perceptivity_v03_scscaled.json  cell_types, modules, A_max, cell_type_class

Method:
  1. Load HPA gene-level nCPM (≈3M rows from rna_single_cell_type.tsv)
  2. Build AnnData (cell types × genes)
  3. sc.pp.normalize_total(target_sum=1e4) + sc.pp.log1p
  4. Score 44 modules with mean_raw (core / activity / etc)
  5. Save as drop-in replacement for the pseudobulk reference

Loaded via: pct.load_hpa_perceptivity(mode='single_cell_scaled')
"""

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from perceptome.catalog import load_catalog, get_genes  # noqa: E402

HPA_DIR = Path("/Users/teo/Desktop/research/perceptual_modules/paper2/data")
HPA_TSV = HPA_DIR / "rna_single_cell_type.tsv"
HPA_CTCLASS = HPA_DIR / "rna_single_cell_type_cell_types.tsv"

OUT_DIR = Path(__file__).resolve().parents[1] / "perceptome" / "perceptivity" / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_NPZ = OUT_DIR / "hpa_perceptivity_v03_scscaled.npz"
OUT_META = OUT_DIR / "hpa_perceptivity_v03_scscaled.json"


def main():
    print(f"Loading catalog (44 modules expected)...")
    cat = load_catalog()
    assert cat["n_modules"] == 44, cat["n_modules"]
    modules = sorted(cat["modules"].keys())

    core_genes = {m: set(get_genes(m, "core", cat)) for m in modules}
    activity_genes = {m: set(get_genes(m, "activity", cat)) for m in modules}

    print(f"Loading HPA per-gene per-cell-type nCPM matrix from {HPA_TSV}...")
    # Streaming parse → dict {(cell_type, gene): nCPM}
    gene_set = set()
    cell_type_set = set()
    rows = []
    with open(HPA_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            ct = r["Cell type"]
            gene = r["Gene name"]
            ncpm = float(r["nCPM"])
            rows.append((ct, gene.upper(), ncpm))
            gene_set.add(gene.upper())
            cell_type_set.add(ct)
    print(f"  parsed {len(rows)} rows, {len(gene_set)} genes, {len(cell_type_set)} cell types")

    cell_types = sorted(cell_type_set)
    genes = sorted(gene_set)
    ct_index = {c: i for i, c in enumerate(cell_types)}
    g_index = {g: i for i, g in enumerate(genes)}

    print(f"Building expression matrix ({len(cell_types)} × {len(genes)})...")
    X = np.zeros((len(cell_types), len(genes)), dtype=np.float32)
    for ct, gene, ncpm in rows:
        X[ct_index[ct], g_index[gene]] = ncpm

    a = ad.AnnData(X=X,
                   obs=pd.DataFrame(index=cell_types),
                   var=pd.DataFrame(index=genes))
    print(f"  AnnData: {a.shape}")

    print("Applying single-cell preprocessing (normalize_total 1e4 + log1p)...")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    print("Scoring 44 modules (core and activity gene sets) per cell type...")
    X_dense = a.X if isinstance(a.X, np.ndarray) else a.X.toarray()
    var_index = {g: i for i, g in enumerate(a.var_names)}

    R = np.zeros((len(cell_types), len(modules)), dtype=np.float32)
    A = np.zeros((len(cell_types), len(modules)), dtype=np.float32)
    for j, m in enumerate(modules):
        # Readiness: mean over core genes present
        present_core = [var_index[g] for g in core_genes[m] if g in var_index]
        if present_core:
            R[:, j] = X_dense[:, present_core].mean(axis=1)
        # Activity: mean over activity-target genes present
        present_act = [var_index[g] for g in activity_genes[m] if g in var_index]
        if present_act:
            A[:, j] = X_dense[:, present_act].mean(axis=1)

    C = R - A
    A_max = A.max(axis=0)
    headroom = A_max[None, :] - A

    print(f"  R range: [{R.min():.3f}, {R.max():.3f}]  (was 0–4 in pseudobulk)")
    print(f"  A range: [{A.min():.3f}, {A.max():.3f}]")
    print(f"  μ per module: range [{R.mean(axis=0).min():.3f}, {R.mean(axis=0).max():.3f}]")
    print(f"  σ per module: range [{R.std(axis=0).min():.3f}, {R.std(axis=0).max():.3f}]")

    # Load cell-type-class mapping for metadata parity with pseudobulk reference
    ct_class = {}
    if HPA_CTCLASS.exists():
        with open(HPA_CTCLASS) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                ct_class[r["Cell type"]] = r.get("Cell type class", "")

    print(f"Saving to {OUT_NPZ} and {OUT_META}...")
    np.savez_compressed(OUT_NPZ, R=R, A=A, C=C, headroom=headroom)
    meta = {
        "cell_types": cell_types,
        "modules": modules,
        "A_max": A_max.tolist(),
        "cell_type_class": [ct_class.get(c, "") for c in cell_types],
        "mode": "single_cell_scaled",
        "preprocessing": "normalize_total(target_sum=1e4) + log1p applied to HPA pseudobulk nCPM",
        "catalog_version": "0.3",
        "version": "0.3.scscaled",
    }
    with open(OUT_META, "w") as f:
        json.dump(meta, f, indent=2)

    print("Done.")


if __name__ == "__main__":
    main()
