"""Recompute aging reference vectors on v0.3 eigenspace (44 modules, 9 PCs).

Inputs (already on disk):
  blood        census_full_blood.h5ad         (2,607 cells, age_bin = old/young)
  bone_marrow  census_full_bone_marrow.h5ad   (3,000 cells, age_bin = old/young)

Procedure (matches v0.1 aging vector construction):
  1. Load h5ad, ensure gene-symbol var_names (CELLxGENE → feature_name)
  2. Score 44 modules (mean_raw, gene_set='core')
  3. Compute mean(old) − mean(young) per module → tissue-specific aging direction
  4. Project into v0.3 9-PC eigenspace
  5. Save to perceptome/reference/data/aging_reference.json

Output schema (compatible with compare/references.py loader):
  inflammaging_direction  blood old−young direction (9-PC vector)
  collapse_direction      bone marrow old−young direction (9-PC vector)
  tissue_vectors          per-tissue per-module deltas + eigenspace projection
  shared_modules          modules where both tissues show same-sign |delta| > 0.05
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scanpy as sc  # noqa: E402

import perceptome as pct  # noqa: E402

CENSUS_DIR = Path("/Users/teo/Desktop/research/paper1/oscilatory/data/census_full")
OUT = Path(__file__).resolve().parents[1] / "perceptome" / "reference" / "data" / "aging_reference.json"

TISSUES = {
    "blood": {
        "h5ad": CENSUS_DIR / "census_full_blood.h5ad",
        "pattern_label": "inflammaging",
    },
    "bone_marrow": {
        "h5ad": CENSUS_DIR / "census_full_bone_marrow.h5ad",
        "pattern_label": "collapse",
    },
}

SHARED_THRESHOLD = 0.05  # min |delta| for "shared" classification


def _ensure_gene_symbols(adata):
    sample = list(adata.var_names[:5])
    if all(s.isdigit() for s in sample):
        if "feature_name" in adata.var.columns:
            adata.var_names = adata.var["feature_name"].astype(str).values
            adata.var_names_make_unique()
            return True
    return False


def _normalize_if_needed(adata):
    X = adata.X
    sample = X[:100].toarray() if hasattr(X, "toarray") else X[:100]
    if sample.max() > 100:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return True
    return False


def process_tissue(name, cfg, eigenspace_loadings, module_order):
    print(f"\n{'='*70}\n{name.upper()}: {cfg['h5ad']}\n{'='*70}")
    if not cfg["h5ad"].exists():
        print(f"  SKIP — file missing")
        return None

    adata = sc.read_h5ad(cfg["h5ad"])
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    if _ensure_gene_symbols(adata):
        print(f"  Remapped var_names to gene symbols")
    if _normalize_if_needed(adata):
        print(f"  Normalized + log1p applied")

    if "age_bin" not in adata.obs.columns:
        print(f"  ERROR: no age_bin column")
        return None

    age_vals = list(adata.obs["age_bin"].unique())
    print(f"  age_bin values: {age_vals}")
    if "old" not in age_vals or "young" not in age_vals:
        print(f"  ERROR: need 'old' and 'young' values")
        return None

    n_old = (adata.obs["age_bin"] == "old").sum()
    n_young = (adata.obs["age_bin"] == "young").sum()
    print(f"  old: {n_old}, young: {n_young}")

    print(f"  Scoring 44 modules (mean_raw, gene_set=core)...")
    scores = pct.score_modules(adata, method="mean_raw", gene_set="core")["scores"]

    old_mask = adata.obs["age_bin"] == "old"
    young_mask = adata.obs["age_bin"] == "young"

    delta = scores.loc[adata.obs.index[old_mask]].mean() - scores.loc[adata.obs.index[young_mask]].mean()
    delta_vec = np.array([delta.get(m, 0.0) for m in module_order])
    eig_vec = delta_vec @ eigenspace_loadings

    top_up = delta.sort_values(ascending=False).head(5)
    top_down = delta.sort_values().head(5)
    print(f"  top +Δ modules: {dict(top_up.round(3))}")
    print(f"  top −Δ modules: {dict(top_down.round(3))}")
    print(f"  eigenspace magnitude: {np.linalg.norm(eig_vec):.4f}")
    print(f"  PC1 component: {eig_vec[0]:+.4f}")

    return {
        "name": name,
        "pattern_label": cfg["pattern_label"],
        "n_old": int(n_old),
        "n_young": int(n_young),
        "module_deltas": {m: float(d) for m, d in zip(module_order, delta_vec)},
        "eigenspace_vector": eig_vec.tolist(),
        "magnitude": float(np.linalg.norm(eig_vec)),
        "pc1_direction": "up" if eig_vec[0] > 0 else "down",
    }


def main():
    print("Recomputing aging reference on v0.3 eigenspace (44 modules, 9 PCs)\n")

    from perceptome.eigenspace.project import _load_reference
    ref = _load_reference()
    module_order = ref["module_order"]
    n_pcs = ref["n_pcs"]
    eigenspace_loadings = np.zeros((len(module_order), n_pcs))
    for j in range(n_pcs):
        for i, m in enumerate(module_order):
            eigenspace_loadings[i, j] = ref["loadings"][f"PC{j+1}"].get(m, 0.0)

    tissue_results = {}
    for name, cfg in TISSUES.items():
        result = process_tissue(name, cfg, eigenspace_loadings, module_order)
        if result is not None:
            tissue_results[name] = result

    if "blood" not in tissue_results or "bone_marrow" not in tissue_results:
        print("\nERROR: need both tissues for shared-module computation")
        sys.exit(1)

    # Shared modules: same sign + |delta| > threshold in both tissues
    blood_d = tissue_results["blood"]["module_deltas"]
    bm_d = tissue_results["bone_marrow"]["module_deltas"]
    shared = []
    for m in module_order:
        b = blood_d[m]
        x = bm_d[m]
        if abs(b) > SHARED_THRESHOLD and abs(x) > SHARED_THRESHOLD and (b > 0) == (x > 0):
            shared.append(m)
    print(f"\nShared aging modules ({len(shared)}, |Δ|>{SHARED_THRESHOLD}, same sign in both):")
    print(f"  {shared}")

    out = {
        "version": "0.3",
        "note": "Recomputed 2026-05-10 on v0.3 eigenspace (44 modules, 9 PCs). Aging is tissue-specific, not universal PC1+: blood shows inflammaging signature (PC1+), bone marrow shows collapse signature (PC1−). Source: CELLxGENE Census full_blood + full_bone_marrow with age_bin annotation.",
        "regenerated": "2026-05-10",
        "eigenspace_n_pcs": 9,
        "n_modules": 44,
        "scoring_method": "mean_raw",
        "gene_set": "core",
        "supersedes": "v0.1 aging_reference.json (12-PC), pinned as aging_reference_v02pinned.json",
        "tissue_vectors": {
            name: {
                "vector_PCs": r["eigenspace_vector"],
                "module_deltas": r["module_deltas"],
                "magnitude": r["magnitude"],
                "pc1_direction": r["pc1_direction"],
                "pattern": r["pattern_label"],
                "n_old": r["n_old"],
                "n_young": r["n_young"],
            }
            for name, r in tissue_results.items()
        },
        "shared_modules": shared,
        "inflammaging_direction": tissue_results["blood"]["eigenspace_vector"],
        "collapse_direction": tissue_results["bone_marrow"]["eigenspace_vector"],
    }

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\n{'='*70}\nwrote {OUT}\n{'='*70}")


if __name__ == "__main__":
    main()
