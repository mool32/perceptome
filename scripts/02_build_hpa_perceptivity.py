"""Precompute HPA perceptivity reference (154 cell types × 44 modules).

Inputs:
  HPA_DIR/rna_single_cell_type.tsv             gene-level nCPM (≈3M rows)
  HPA_DIR/rna_single_cell_type_cell_types.tsv  cell type → "Cell type class"

Outputs (under perceptome/perceptivity/data/):
  hpa_perceptivity_v03.npz   R, A, C, headroom matrices (154 × 44, float32)
  hpa_perceptivity_v03.json  cell_types, modules, A_max, cell_type_class

Method:
  R(C, M) = mean_{g in core_genes(M) found in HPA} log1p(nCPM(g, C))
  A(C, M) = mean_{g in activity_genes(M) found in HPA} log1p(nCPM(g, C))
  C       = R - A
  headroom = A_max(M) − A(C, M)

This matches v1.4 metric design (continuous, log1p(nCPM), per-cell-type mean).
"""

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

# Make the in-tree package importable when run as a script
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from perceptome.catalog import load_catalog, get_genes  # noqa: E402

HPA_DIR = Path("/Users/teo/Desktop/research/paper2/data")
HPA_TSV = HPA_DIR / "rna_single_cell_type.tsv"
HPA_CTCLASS = HPA_DIR / "rna_single_cell_type_cell_types.tsv"

OUT_DIR = Path(__file__).resolve().parents[1] / "perceptome" / "perceptivity" / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_NPZ = OUT_DIR / "hpa_perceptivity_v03.npz"
OUT_META = OUT_DIR / "hpa_perceptivity_v03.json"


def main():
    print(f"Loading catalog (44 modules expected)...")
    cat = load_catalog()
    assert cat["n_modules"] == 44, cat["n_modules"]
    modules = sorted(cat["modules"].keys())

    core_genes = {m: set(get_genes(m, "core", cat)) for m in modules}
    activity_genes = {m: set(get_genes(m, "activity", cat)) for m in modules}

    all_genes = set().union(*core_genes.values(), *activity_genes.values())
    print(f"  {len(modules)} modules, {len(all_genes)} unique genes referenced")

    print(f"Loading HPA cell type class map...")
    ct_class = {}
    with open(HPA_CTCLASS) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ct_class[row["Cell type"]] = row["Cell type class"]
    print(f"  {len(ct_class)} cell types")

    print(f"Streaming HPA expression file: {HPA_TSV}")
    print(f"  (3M rows; only retaining ~{len(all_genes)} genes used by catalog)")

    expr = defaultdict(dict)  # gene_name -> {cell_type: nCPM}
    n_rows = 0
    n_kept = 0
    with open(HPA_TSV) as f:
        next(f)  # header
        for line in f:
            n_rows += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            gname = parts[1]
            if gname not in all_genes:
                continue
            ct = parts[2]
            try:
                ncpm = float(parts[3])
            except ValueError:
                continue
            expr[gname][ct] = ncpm
            n_kept += 1
    print(f"  scanned {n_rows} rows, kept {n_kept}; {len(expr)} catalog genes have HPA values")

    cell_types = sorted({ct for d in expr.values() for ct in d})
    print(f"  {len(cell_types)} HPA cell types in expression data")

    def gene_set_mean_log(genes, ct):
        vals = [math.log1p(expr[g][ct]) for g in genes if g in expr and ct in expr[g]]
        return float(np.mean(vals)) if vals else float("nan")

    n_ct = len(cell_types)
    n_mod = len(modules)
    R = np.full((n_ct, n_mod), np.nan, dtype=np.float32)
    A = np.full((n_ct, n_mod), np.nan, dtype=np.float32)

    print(f"Computing R (core) and A (activity) for {n_ct} × {n_mod} = {n_ct * n_mod} cells...")
    for j, m in enumerate(modules):
        cg = core_genes[m]
        ag = activity_genes[m]
        for i, ct in enumerate(cell_types):
            R[i, j] = gene_set_mean_log(cg, ct)
            A[i, j] = gene_set_mean_log(ag, ct)

    C = R - A
    A_max = np.nanmax(A, axis=0)
    headroom = A_max[None, :] - A

    np.savez_compressed(OUT_NPZ, R=R, A=A, C=C, headroom=headroom)
    print(f"  wrote {OUT_NPZ} ({OUT_NPZ.stat().st_size / 1024:.1f} KB)")

    meta = {
        "version": "0.3",
        "n_cell_types": n_ct,
        "n_modules": n_mod,
        "cell_types": cell_types,
        "modules": modules,
        "A_max": A_max.tolist(),
        "cell_type_class": [ct_class.get(ct, "unknown") for ct in cell_types],
        "method": "mean log1p(nCPM)",
        "source_hpa_file": str(HPA_TSV),
        "catalog_version": cat["version"],
    }
    with open(OUT_META, "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print(f"  wrote {OUT_META}")

    # Sanity: a few known reference cell types
    print("\nSanity: NPAS4 (neuron module) per cell-type-class group:")
    j_npas4 = modules.index("NPAS4")
    by_class = defaultdict(list)
    for i, ct in enumerate(cell_types):
        by_class[ct_class.get(ct, "unknown")].append(A[i, j_npas4])
    for klass in sorted(by_class):
        vals = [v for v in by_class[klass] if not np.isnan(v)]
        if vals:
            print(f"  {klass:30s} mean A_NPAS4 = {np.mean(vals):.2f} (n={len(vals)})")

    print("\nSanity: UPR-ATF6 (saturation reference) per cell-type-class group:")
    j_atf6 = modules.index("UPR-ATF6")
    by_class = defaultdict(list)
    for i, ct in enumerate(cell_types):
        by_class[ct_class.get(ct, "unknown")].append(A[i, j_atf6])
    for klass in sorted(by_class):
        vals = [v for v in by_class[klass] if not np.isnan(v)]
        if vals:
            print(f"  {klass:30s} mean A_UPR-ATF6 = {np.mean(vals):.2f} (n={len(vals)})")


if __name__ == "__main__":
    main()
