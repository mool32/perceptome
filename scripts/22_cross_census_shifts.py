#!/usr/bin/env python3
"""Phase 2b/#1 — extend the signed tumor−normal shifts using CELLxGENE-Census
normal cell-of-origin atlases (cross-accession, same modality/pipeline as the
malignant caches → HPA calibration cancels in the difference; residual confound =
cross-accession batch + subsample). TIER 2 (approximate), distinct from the
within-dataset gold-standard (TIER 1).

Only tissues whose parenchymal cell-of-origin is actually present in the (immune-
biased ~3000-cell) census_full subsample qualify.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
ROOT = Path("/Users/teo/Desktop/research/perceptual_modules")
TOOL = ROOT / "tool" / "perceptome2"
CEN = ROOT / "paper1/oscilatory/data/census_full"
sys.path.insert(0, str(TOOL))
import perceptome as pct  # noqa: E402
from perceptome.eigenspace.project import _load_reference  # noqa: E402

ref = _load_reference(None); mo = ref["module_order"]; npc = ref["n_pcs"]
std_ref = np.array(ref["std_per_module"]); std_ref[std_ref < 1e-10] = 1.0
loadings = np.array([[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(npc)] for m in mo])
cc = json.load(open(ROOT / "transition_vectors/results/per_pc_decomposition.json"))["class_centroids"]
Mtr = np.array([cc[n] for n in cc]); Mtr = Mtr / np.linalg.norm(Mtr, axis=1, keepdims=True)
ms = np.linalg.svd(Mtr, full_matrices=False)[2][0]
if (Mtr[list(cc).index("TUMOR_TME")] @ ms) < 0:
    ms = -ms
_core = set(pct.get_genes("HSF1", "core") + pct.get_genes("p53", "core"))

def cos(a, b): return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))

def centroid(adata, mask=None):
    if "feature_name" in adata.var.columns and len(_core & set(adata.var_names)) == 0:
        adata.var_names = adata.var["feature_name"].astype(str).values; adata.var_names_make_unique()
    if mask is not None:
        adata = adata[mask]
    if float(adata.X[:200].max()) > 50:
        sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    return pct.score_modules(adata)["scores"].reindex(columns=mo).mean().values, adata.n_obs

# TIER 1 — within-dataset gold-standard (from script 21)
GOLD = {"HCC": (+0.47, 0.746), "BRCA": (-0.62, 0.718), "BCC": (+0.85, 0.254)}

# TIER 2 — cross-Census cell-of-origin
CROSS = {
    "RCC": ("results/cancer_scaled/cache/RCC_clearcell_malignant.h5ad", "kidney", "epithelial cell of proximal tubule"),
    "PDAC": ("results/cancer_v6_5_extras/PDAC_Peng2019_malignant.h5ad", "pancreas", "pancreatic ductal cell"),
    "Glioblastoma": ("results/cancer_scaled/cache/Glioblastoma_malignant.h5ad", "brain", "radial glial cell"),
    "Neuroblastoma": ("results/cancer_v6_4_extras/Neuroblastoma_malignant.h5ad", "brain", "neuroblast (sensu Vertebrata)"),
    "DLBCL": ("results/cancer_v6_4_extras/DLBCL_malignant.h5ad", "spleen", "B cell"),
}
cross = {}
for name, (mp, tissue, origin) in CROSS.items():
    mc, nm = centroid(sc.read_h5ad(ROOT / mp))
    atlas = sc.read_h5ad(CEN / f"census_full_{tissue}.h5ad")
    nmask = atlas.obs["cell_type"].astype(str).values == origin
    nc, nn = centroid(atlas, nmask)
    shift = mc - nc
    cross[name] = (cos((shift / std_ref) @ loadings, ms), float(np.linalg.norm(shift)), nm, nn, origin)
    print(f"  {name:14s} malig n={nm:6d}  normal '{origin}' n={nn:5d}  cos={cross[name][0]:+.2f}  ‖{cross[name][1]:.2f}‖")

# ---- combined report ----
print("\nSigned tumor−normal shift on the main-sequence axis (cos>0 = chronic-ward):")
print("  TIER 1 (within-dataset, calibration-clean):")
for k, (c, m) in GOLD.items():
    print(f"    {k:14s} cos={c:+.2f}  ‖{m:.2f}‖")
print("  TIER 2 (cross-Census cell-of-origin, approximate — batch+subsample confound):")
for k, (c, m, *_ ) in cross.items():
    print(f"    {k:14s} cos={c:+.2f}  ‖{m:.2f}‖")
allc = {**{k: v[0] for k, v in GOLD.items()}, **{k: v[0] for k, v in cross.items()}}
print(f"\n  {sum(v > 0 for v in allc.values())}/{len(allc)} chronic-ward, "
      f"{sum(v < 0 for v in allc.values())}/{len(allc)} acute-ward → direction is SPLIT (tissue-dependent).")

# ---- figure ----
fig, ax = plt.subplots(figsize=(9, 6))
items = sorted(allc.items(), key=lambda kv: kv[1])
for y, (name, c) in enumerate(items):
    tier1 = name in GOLD
    col = "crimson" if c > 0 else "steelblue"
    ax.barh(y, c, color=col, alpha=0.85 if tier1 else 0.4,
            edgecolor="k", linewidth=0.8 if tier1 else 0.3,
            hatch=None if tier1 else "///")
    ax.text(c + (0.02 if c > 0 else -0.02), y, name, va="center",
            ha="left" if c > 0 else "right", fontsize=9)
ax.axvline(0, color="k", lw=1)
ax.set_yticks([])
ax.set_xlim(-0.9, 1.0)
ax.set_xlabel("cos(tumor−normal shift, main-sequence axis)\n◀ acute-ward            chronic-ward ▶")
ax.set_title("Tumor−normal shift direction is tissue-dependent\n"
             "(solid = within-dataset gold-standard; hatched = cross-Census, approximate)")
ax.scatter([], [], marker="s", c="crimson", label="chronic-ward"); ax.scatter([], [], marker="s", c="steelblue", label="acute-ward")
ax.legend(loc="lower right", fontsize=9)
plt.tight_layout()
out = TOOL / "scripts" / "fig_shift_directions.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nfigure → {out}")
