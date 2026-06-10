#!/usr/bin/env python3
"""Phase 2b/#1 — fully non-circular cancer trajectory: recompute the Sun2021 HCC
paired tumor−normal shift from raw scRNA and project it onto the main-sequence axis.

Non-circular because: (a) computed from real paired tumor/normal cells (not from
the attractor-cell definition), (b) projected in MODULE space normalized by the
reference per-module std, so it never touches the cancer-derived attractor.
Confirms (or not) that cancer sits at the chronic pole of the main sequence,
independently of the transition_vectors TUMOR_TME class (+0.85).
"""
import json
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

ROOT = Path("/Users/teo/Desktop/research/perceptual_modules")
TOOL = ROOT / "tool" / "perceptome2"
sys.path.insert(0, str(TOOL))
import perceptome as pct  # noqa: E402
from perceptome.eigenspace.project import _load_reference  # noqa: E402

ref = _load_reference(None)
mo = ref["module_order"]; npc = ref["n_pcs"]
std_ref = np.array(ref["std_per_module"]); std_ref[std_ref < 1e-10] = 1.0
loadings = np.array([[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(npc)] for m in mo])

# ---- main-sequence axis (#1) from transition-class centroids; R1 (#2) from varimax ----
cc = json.load(open(ROOT / "transition_vectors/results/per_pc_decomposition.json"))["class_centroids"]
Mtr = np.array([cc[n] for n in cc]); Mtr = Mtr / np.linalg.norm(Mtr, axis=1, keepdims=True)
ms = np.linalg.svd(Mtr, full_matrices=False)[2][0]
if (Mtr @ ms).sum() < 0:
    ms = -ms
# orient ms so the cancer transition classes are POSITIVE (chronic pole = +)
if (Mtr[list(cc).index("TUMOR_TME")] @ ms) < 0:
    ms = -ms
rot6 = pct.varimax(loadings[:, :6])[0]
R1_9 = np.zeros(npc); R1_9[:6] = pct.varimax(loadings[:, :6])[1][:, 0]

def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))

# ---- load + score HCC ----
print("loading Sun2021 HCC h5ad (71915 × 25712)...")
a = sc.read_h5ad(ROOT / "results/cancer_sun2021_hcc/sun2021.h5ad")
if float(a.X[:200].max()) > 50:
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    print("  applied normalize_total(1e4)+log1p (raw counts detected)")
scores = pct.score_modules(a)["scores"].reindex(columns=mo)

def shift_vec(cell_type):
    m = a.obs["celltype"].values == cell_type
    tum = m & (a.obs["site"].values == "Tumor")
    nor = m & (a.obs["site"].values == "Normal")
    if tum.sum() < 20 or nor.sum() < 20:
        return None, tum.sum(), nor.sum()
    d = scores[tum].mean().values - scores[nor].mean().values     # 44-module shift
    return (d / std_ref) @ loadings, tum.sum(), nor.sum()          # → 9-PC

print("\nTumor−normal shift projected onto the main-sequence axis (chronic = +):")
print(f"  {'compartment':14s} {'n_tum':>6s} {'n_norm':>6s}  cos(shift, main-seq)  cos(shift, R1)")
for ct in ["Hepatocyte", "Fibroblast", "Myeloid", "Endothelial", "T/NK", "B"]:
    v, nt, nn = shift_vec(ct)
    if v is None:
        print(f"  {ct:14s} {nt:>6d} {nn:>6d}  (too few cells)"); continue
    tag = "  ← MALIGNANT" if ct == "Hepatocyte" else ""
    print(f"  {ct:14s} {nt:>6d} {nn:>6d}        {cos(v, ms):+.2f}             {cos(v, R1_9):+.2f}{tag}")

print(f"\nreference: transition_vectors TUMOR_TME sits at +0.85 on this axis; "
      f"a positive cos here independently confirms cancer → chronic pole.")
