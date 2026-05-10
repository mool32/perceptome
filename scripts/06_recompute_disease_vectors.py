"""Recompute disease perturbation vectors on v0.3 eigenspace (44 modules, 9 PCs).

Inputs (already on disk under perceptual_modules/results/diseases/):
  D1_RA   ra_cellxgene.h5ad
  D2_AD   GSE138852_grubman2019_AD_entorhinal.h5ad
  D3_IPF  ipf_habermann2020.h5ad
  D4_DKD  dkd_cellxgene.h5ad

Procedure (matches v0.1 disease pipeline, only catalog + eigenspace updated):
  1. Load h5ad, normalize+log1p if needed
  2. For each cell-type group: compute disease-cell scores − control-cell scores
     via score_modules(method='mean_raw', gene_set='core')
  3. Project per-group delta into v0.3 9-PC eigenspace
  4. Save as disease_vectors.json (compatible with compare/references.py loader)

Output:
  perceptome/reference/data/disease_vectors.json
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scanpy as sc  # noqa: E402

import perceptome as pct  # noqa: E402

DISEASES_DIR = Path("/Users/teo/Desktop/research/perceptual_modules/results/diseases")
OUT = Path(__file__).resolve().parents[1] / "perceptome" / "reference" / "data" / "disease_vectors.json"

# ── Cell-group definitions per disease (locked from v0.1 disease_D[1-4]_*.py) ──

D1_RA_GROUPS = {
    "monocytes": ["classical monocyte", "non-classical monocyte"],
    "CD4_T": ["central memory CD4-positive, alpha-beta T cell",
              "naive thymus-derived CD4-positive, alpha-beta T cell",
              "effector memory CD4-positive, alpha-beta T cell",
              "CD4-positive, alpha-beta T cell"],
    "CD8_T": ["naive thymus-derived CD8-positive, alpha-beta T cell",
              "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
              "CD8-positive, alpha-beta memory T cell"],
    "NK": ["natural killer cell"],
    "B_cells": ["naive B cell", "memory B cell"],
    "dendritic": ["myeloid dendritic cell"],
}
D2_AD_GROUPS = {
    "neurons": ["neuron"],
    "microglia": ["mg"],
    "astrocytes": ["astro"],
    "oligodendrocytes": ["oligo"],
    "OPCs": ["OPC"],
}
D3_IPF_GROUPS = {
    "fibroblasts": ["PLIN2+ Fibroblasts", "HAS1 High Fibroblasts", "Myofibroblasts", "Fibroblasts"],
    "macrophages": ["Macrophages"],
    "AT2_epithelial": ["AT2", "Transitional AT2"],
    "basal_epithelial": ["Basal"],
    "endothelial": ["Endothelial Cells"],
    "T_cells": ["T Cells"],
    "monocytes": ["Monocytes"],
    "ciliated": ["Ciliated", "Differentiating Ciliated"],
    "NK_cells": ["NK Cells"],
}
D4_DKD_GROUPS = {
    "proximal_tubular": ["epithelial cell of proximal tubule"],
    "distal_tubular": ["kidney distal convoluted tubule epithelial cell"],
    "thick_ascending_limb": ["kidney loop of Henle thick ascending limb epithelial cell"],
    "collecting_duct": ["renal principal cell", "renal alpha-intercalated cell",
                        "renal beta-intercalated cell"],
    "podocytes": ["podocyte"],
    "endothelial": ["endothelial cell"],
    "parietal_epithelial": ["parietal epithelial cell"],
    "mesangial": ["mesangial cell"],
    "immune": ["leukocyte"],
}

DISEASES = {
    "RA": {
        "h5ad": DISEASES_DIR / "D1_RA" / "ra_cellxgene.h5ad",
        "groups": D1_RA_GROUPS,
        "disease_col": "disease",
        "disease_label": ("rheumatoid arthritis", "Rheumatoid arthritis", "RA"),
        "control_label": ("normal", "Normal", "healthy", "control"),
        "celltype_col": "cell_type",
        "dataset_name": "CELLxGENE RA PBMC (Zhang et al.)",
    },
    "AD": {
        "h5ad": DISEASES_DIR / "D2_AD" / "GSE138852_grubman2019_AD_entorhinal.h5ad",
        "groups": D2_AD_GROUPS,
        "disease_col": None,        # AD uses oupSample-id-encoded labels
        "disease_label": "AD",
        "control_label": "ct",
        "celltype_col": None,        # auto-detect
        "dataset_name": "Grubman 2019 GSE138852 AD entorhinal",
    },
    "IPF": {
        "h5ad": DISEASES_DIR / "D3_IPF" / "ipf_habermann2020.h5ad",
        "groups": D3_IPF_GROUPS,
        "disease_col": "Diagnosis",  # Habermann uses 'Diagnosis' column
        "disease_label": ("IPF", "ipf"),
        "control_label": ("Control", "control", "Donor", "donor"),
        "celltype_col": "celltype",
        "dataset_name": "Habermann 2020 IPF",
    },
    "DKD": {
        "h5ad": DISEASES_DIR / "D4_DKD" / "dkd_cellxgene.h5ad",
        "groups": D4_DKD_GROUPS,
        "disease_col": None,
        "disease_label": ("type 2 diabetes mellitus", "diabetic kidney disease", "DKD"),
        "control_label": ("normal", "Normal", "healthy"),
        "celltype_col": "cell_type",
        "dataset_name": "CELLxGENE DKD",
    },
}


def _auto_detect_col(adata, candidates):
    for c in candidates:
        if c in adata.obs.columns:
            return c
    return None


def _resolve_label(values, candidates):
    for c in candidates if isinstance(candidates, (list, tuple)) else [candidates]:
        if c in values:
            return c
    # fuzzy
    for c in candidates if isinstance(candidates, (list, tuple)) else [candidates]:
        for v in values:
            if c.lower() in str(v).lower():
                return v
    return None


def _normalize_if_needed(adata):
    X = adata.X
    if hasattr(X, "toarray"):
        sample = X[:100].toarray()
    else:
        sample = X[:100]
    if sample.max() > 100:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return True
    return False


def _ensure_gene_symbols(adata):
    """CELLxGENE Census stores var_names as integer indices and gene symbols in
    var['feature_name']. Detect and rewrite var_names so module scoring matches.
    """
    sample = list(adata.var_names[:5])
    if all(s.isdigit() for s in sample):
        if "feature_name" in adata.var.columns:
            adata.var_names = adata.var["feature_name"].astype(str).values
            adata.var_names_make_unique()
            return True
    return False


def process_disease(name, cfg):
    print(f"\n{'='*70}\n{name}: {cfg['h5ad']}\n{'='*70}")

    if not cfg["h5ad"].exists():
        print(f"  SKIP — file missing")
        return None

    adata = sc.read_h5ad(cfg["h5ad"])
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # CELLxGENE Census: integer var_names + gene symbols in var['feature_name']
    if _ensure_gene_symbols(adata):
        print(f"  Remapped var_names: integer→feature_name (gene symbols)")

    # Normalize if needed
    normed = _normalize_if_needed(adata)
    print(f"  {'Normalized + log1p applied' if normed else 'Already log-normalized'}")

    # Resolve cell type column
    ct_col = cfg["celltype_col"] or _auto_detect_col(
        adata, ["cell_type", "celltype", "Celltype", "cellType", "Cell_Type"]
    )
    if ct_col is None:
        print(f"  ERROR: no cell type column found; obs columns: {list(adata.obs.columns)[:10]}")
        return None
    print(f"  cell type column: {ct_col}")

    # Resolve disease/control column
    disease_col = cfg["disease_col"] or _auto_detect_col(
        adata, ["disease", "Disease", "condition", "Condition", "diagnosis", "group", "sample_type"]
    )
    if disease_col is None:
        print(f"  ERROR: no disease column found; obs columns: {list(adata.obs.columns)[:10]}")
        return None
    print(f"  disease column: {disease_col}")

    disease_vals = adata.obs[disease_col].unique()
    print(f"  unique values: {list(disease_vals)[:5]}{'...' if len(disease_vals) > 5 else ''}")

    disease_lbl = _resolve_label(set(disease_vals), cfg["disease_label"])
    control_lbl = _resolve_label(set(disease_vals), cfg["control_label"])
    if disease_lbl is None or control_lbl is None:
        print(f"  ERROR: could not resolve labels — disease={disease_lbl}, control={control_lbl}")
        return None
    print(f"  disease label: '{disease_lbl}' | control label: '{control_lbl}'")

    # Score 44 modules (mean_raw, gene_set='core' for readiness)
    print(f"  Scoring 44 modules (mean_raw, gene_set=core)...")
    scores = pct.score_modules(adata, method="mean_raw", gene_set="core")["scores"]

    per_group = {}
    eigenspace_loadings = None  # cache after first build
    n_pcs = None

    from perceptome.eigenspace.project import _load_reference
    ref = _load_reference()
    module_order = ref["module_order"]
    n_pcs = ref["n_pcs"]
    eigenspace_loadings = np.zeros((len(module_order), n_pcs))
    for j in range(n_pcs):
        for i, m in enumerate(module_order):
            eigenspace_loadings[i, j] = ref["loadings"][f"PC{j+1}"].get(m, 0.0)

    aggregate_delta = np.zeros(len(module_order))
    n_groups_with_data = 0

    for group_name, cell_types in cfg["groups"].items():
        mask = adata.obs[ct_col].isin(cell_types)
        if mask.sum() < 50:
            print(f"    {group_name}: only {mask.sum()} cells — skip")
            continue

        d_mask = mask & (adata.obs[disease_col] == disease_lbl)
        c_mask = mask & (adata.obs[disease_col] == control_lbl)
        if d_mask.sum() < 20 or c_mask.sum() < 20:
            print(f"    {group_name}: too few d={d_mask.sum()}/c={c_mask.sum()} — skip")
            continue

        # Reindex via boolean indexing on the score DataFrame
        d_scores = scores.loc[adata.obs.index[d_mask]]
        c_scores = scores.loc[adata.obs.index[c_mask]]
        delta = d_scores.mean() - c_scores.mean()

        delta_vec = np.array([delta.get(m, 0.0) for m in module_order])
        eig_vec = delta_vec @ eigenspace_loadings

        per_group[group_name] = {
            "n_disease": int(d_mask.sum()),
            "n_control": int(c_mask.sum()),
            "module_deltas": {m: float(delta.get(m, 0.0)) for m in module_order},
            "eigenspace_vector": eig_vec.tolist(),
            "magnitude": float(np.linalg.norm(eig_vec)),
        }
        aggregate_delta += delta_vec
        n_groups_with_data += 1
        print(f"    {group_name}: d={d_mask.sum()}, c={c_mask.sum()}, magnitude={per_group[group_name]['magnitude']:.4f}")

    if n_groups_with_data == 0:
        print(f"  ERROR: no groups had enough data for {name}")
        return None

    aggregate_delta /= n_groups_with_data
    aggregate_eig = aggregate_delta @ eigenspace_loadings

    return {
        "dataset": cfg["dataset_name"],
        "n_cells_total": int(adata.n_obs),
        "disease_label": disease_lbl,
        "control_label": control_lbl,
        "n_cell_groups": n_groups_with_data,
        "per_cell_type": per_group,
        "module_delta_mean_across_groups": {m: float(d) for m, d in zip(module_order, aggregate_delta)},
        "mean_vector_12D": aggregate_eig.tolist(),  # name kept for compat with compare/references loader
        "eigenspace_PC_count": n_pcs,
    }


def main():
    print("Recomputing disease vectors on v0.3 eigenspace (44 modules, 9 PCs)\n")
    out = {}
    for name, cfg in DISEASES.items():
        result = process_disease(name, cfg)
        if result is not None:
            out[name] = result

    if not out:
        print("\nERROR: no diseases successfully processed")
        sys.exit(1)

    out["_meta"] = {
        "version": "0.3",
        "eigenspace_n_pcs": 9,
        "n_modules": 44,
        "scoring_method": "mean_raw",
        "gene_set": "core",
        "regenerated": "2026-05-10",
        "source_pipeline": "tool/perceptome2/scripts/06_recompute_disease_vectors.py",
        "supersedes": "v0.1 disease_vectors.json (12-PC), pinned as disease_vectors_v02pinned.json",
    }

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\n{'='*70}\nwrote {OUT}\nDiseases processed: {[k for k in out if not k.startswith('_')]}\n{'='*70}")


if __name__ == "__main__":
    main()
