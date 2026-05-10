"""Build examples/recipe_cancer_attractor.ipynb — focused use-case (5-7 cells)."""

import sys
from pathlib import Path
import nbformat as nbf

OUT = Path(__file__).resolve().parents[1] / "examples" / "recipe_cancer_attractor.ipynb"


def md(text):
    return nbf.v4.new_markdown_cell(text)


def code(src):
    return nbf.v4.new_code_cell(src)


nb = nbf.v4.new_notebook()
nb.cells = [
    md("""# Recipe — cancer attractor analysis

**One question:** does the (tumor − normal) shift in my paired-cohort scRNA-seq align with the cancer-transformation direction identified in Paper 4.2?

**Inputs:** an AnnData with paired tumor and normal cells from the same patient, log-normalized.

**Output:** cosine alignment + PASS/FAIL against the locked Paper 4.2 P3 threshold (>+0.20).

For the full tutorial walkthrough, see [`tutorial.ipynb`](tutorial.ipynb). This recipe is the minimal viable analysis."""),

    code("""import scanpy as sc
import perceptome as pct

# Replace with your paired tumor + adjacent normal scRNA-seq
# adata = sc.read_h5ad("paired_cohort.h5ad")
# adata.obs must have a column distinguishing tumor cells from normal cells

# For demonstration: synthetic split of PBMC3K (NOT real biology)
adata = sc.datasets.pbmc3k_processed()
adata = adata[adata.obs["louvain"].isin(["B cells", "CD4 T cells", "CD8 T cells", "NK cells", "CD14+ Monocytes", "Dendritic cells"])].copy()
adata.obs["arm"] = ["tumor" if i % 2 == 0 else "normal" for i in range(adata.n_obs)]
print(f"{adata.n_obs} cells, arms: {dict(adata.obs.arm.value_counts())}")"""),

    code("""# Score per-cell-type tumor vs normal delta in module space
# Replace 'arm' with your tumor/normal column, 'louvain' with your cell-type column

scores = pct.score_modules(adata, method="mean_raw")["scores"]

# Per-cell-type shift = mean(tumor) − mean(normal)
import pandas as pd
result_rows = []
for ct in adata.obs["louvain"].unique():
    ct_mask = adata.obs["louvain"] == ct
    tumor = scores.loc[adata.obs.index[ct_mask & (adata.obs["arm"] == "tumor")]]
    normal = scores.loc[adata.obs.index[ct_mask & (adata.obs["arm"] == "normal")]]
    if len(tumor) < 10 or len(normal) < 10:
        continue
    shift = (tumor.mean() - normal.mean())  # 44-module Δ vector
    align = pct.attractor_alignment(shift, mode="modules")
    result_rows.append({
        "cell_type": ct,
        "n_tumor": len(tumor),
        "n_normal": len(normal),
        "cosine_with_attractor": align["cosine"],
        "passes_p3_threshold": align["passes_p3_threshold"],
    })
result = pd.DataFrame(result_rows).sort_values("cosine_with_attractor", ascending=False)
print(result.to_string(index=False))"""),

    md("""**Reading the result:**

- `cosine_with_attractor > +0.20` (Paper 4.2 P3 PASS threshold) → cell type's shift aligns with the cancer-transformation direction
- `4 of 6` (or more) cell types passing = cohort-level PASS, matching Sun 2021 paired HCC outcome
- Negative cosines indicate the shift goes *against* the attractor direction — biologically meaningful (some compartments don't transform, e.g. cytotoxic immune cells)

In the synthetic PBMC3K split above the cosines will be near zero (random split, no real biology). On real paired tumor/normal data (e.g. Sun 2021 GSE149614 HCC) you'd expect the malignant compartment + 2-3 stromal compartments to pass."""),

    md("""## Eigenspace alignment (optional — coordinate-space view)

If you want the alignment in eigenspace coordinates instead of module space, project first then use `mode='eigenspace'`:"""),

    code("""# Eigenspace mode: project the per-cell-type shift, then test alignment in 9-PC space
# Useful if you're combining with PC-space visualizations

import numpy as np
shift_vector = result.copy()
for _, row in result.iterrows():
    ct = row["cell_type"]
    ct_mask = adata.obs["louvain"] == ct
    tumor_scores = scores.loc[adata.obs.index[ct_mask & (adata.obs["arm"] == "tumor")]]
    normal_scores = scores.loc[adata.obs.index[ct_mask & (adata.obs["arm"] == "normal")]]
    if len(tumor_scores) < 10 or len(normal_scores) < 10:
        continue
    shift_module = (tumor_scores.mean() - normal_scores.mean()).to_frame().T
    shift_eig = pct.project(shift_module)["coordinates"].iloc[0].values
    align_eig = pct.attractor_alignment(shift_eig, mode="eigenspace")
    print(f"{ct:25s}  cos_eigenspace = {align_eig['cosine']:+.3f}  passes={align_eig['passes_p3_threshold']}")"""),

    md("""## Reference

- **Paper 4.2** identified the 8-cell attractor cluster (gastric chief, pancreatic acinar, parietal, megakaryocytes, cytotrophoblasts, migrating cytotrophoblasts, late primary spermatocytes, gastric progenitor cells) as the convergent destination of 11 cancer types from 11 organ systems.
- **P3 PASS** verdict was: ≥ 4 of 6 cell types in Sun 2021 GSE149614 paired hepatocellular carcinoma cohort show cosine alignment > +0.20.
- **Locked finding:** "Cancer cells move in the attractor direction in capacity space without arriving at the attractor capacity profile" — direction is replicable across transformation, destination correspondence is not.

See `pct.load_attractor_direction()` for the full reference data + provenance."""),
]

nbf.write(nb, OUT)
print(f"wrote {OUT}")
print(f"  cells: {len(nb.cells)}")
