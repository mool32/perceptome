#!/usr/bin/env python3
"""Phase 2b/#1 — the cellular 'main sequence': do distinct biological processes
live as trajectories on ONE shared low-dimensional axis of the manifold?

Assembles every available process-trajectory as a 9-PC direction in the SAME
eigenspace frame (raw module-delta @ loadings, matching transition_vectors):
  * 11 transition-class centroids (transition_vectors/, 81 vectors, pre-registered,
    independent of each other — non-circular)
  * aging: inflammaging + collapse (bundled aging_reference)
  * differentiation-completion: embryonic→terminal HPA gradient (derived here)

Then quantifies: effective dimensionality of the trajectory set, the dominant
shared axis (the 'main sequence'), each process's position on it, and the
module interpretation of that axis.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/Users/teo/Desktop/research/perceptual_modules")
TOOL = ROOT / "tool" / "perceptome2"
sys.path.insert(0, str(TOOL))
import perceptome as pct  # noqa: E402
from perceptome.eigenspace.project import _load_reference  # noqa: E402

ref = _load_reference(None)
mo = ref["module_order"]; npc = ref["n_pcs"]
loadings = np.array([[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(npc)] for m in mo])

def unit(v):
    v = np.asarray(v, float)
    return v / (np.linalg.norm(v) + 1e-12)

# ---- assemble process-trajectory vectors (all 9-PC) ----
vecs, labels, src = [], [], []

cc = json.load(open(ROOT / "transition_vectors/results/per_pc_decomposition.json"))["class_centroids"]
for name, v in cc.items():
    vecs.append(unit(v)); labels.append(name); src.append("transition")

ag = json.load(open(TOOL / "perceptome/reference/data/aging_reference.json"))
vecs.append(unit(ag["inflammaging_direction"])); labels.append("AGING_INFLAMMAGING"); src.append("aging")
vecs.append(unit(ag["collapse_direction"]));     labels.append("AGING_COLLAPSE");     src.append("aging")

# differentiation-completion: embryonic→terminal HPA gradient (raw readiness @ loadings)
R = pct.load_hpa_perceptivity(mode="pseudobulk")["R"].reindex(columns=mo)
cells = list(R.index)
embryonic = ["cytotrophoblasts", "migrating cytotrophoblasts", "syncytiotrophoblasts",
             "extravillous trophoblasts", "oocytes", "early primary spermatocytes",
             "late primary spermatocytes", "early spermatids", "late spermatids",
             "differentiating spermatogonia", "undifferentiated spermatogonia"]
terminal = ["erythrocytes", "cardiomyocytes", "hepatocytes", "brain excitatory neurons",
            "enterocytes", "adipocytes", "rod photoreceptor cells"]
emb_R = R.loc[[c for c in embryonic if c in cells]].mean(0).values
term_R = R.loc[[c for c in terminal if c in cells]].mean(0).values
diff_axis = (term_R - emb_R) @ loadings        # differentiation COMPLETION (toward terminal)
vecs.append(unit(diff_axis)); labels.append("DIFFERENTIATION_COMPLETION"); src.append("derived")

M = np.array(vecs)                              # P × 9, unit rows
P = len(labels)

# ---- effective dimensionality of the trajectory set ----
U, S, Vt = np.linalg.svd(M, full_matrices=False)
var = S ** 2 / (S ** 2).sum()
pr = (S.sum() ** 2) / (S ** 2).sum()
print("=" * 76)
print(f"TRAJECTORY SET: {P} process-vectors in 9-PC space")
print("=" * 76)
print(f"  shared-component variance: {np.round(var[:5], 3)}")
print(f"  participation ratio (effective # shared dims): {pr:.2f}")
print(f"  variance on 1st shared axis (the 'main sequence'): {var[0]:.0%}; on top 2: {var[:2].sum():.0%}")

# ---- the main-sequence axis (1st shared direction) + each process's position ----
ms = Vt[0]
# orient so most processes are positive (sign is arbitrary)
if (M @ ms).sum() < 0:
    ms = -ms
coord = M @ ms
order = np.argsort(coord)
print("\nMAIN-SEQUENCE ORDERING (projection of each process onto the dominant shared axis):")
for i in order:
    print(f"  {coord[i]:+.2f}  {labels[i]:30s} [{src[i]}]")

# ---- module interpretation of the main-sequence axis ----
ms_modules = loadings @ ms          # 9-PC axis → 44 modules
mods = sorted(zip(mo, ms_modules), key=lambda x: -abs(x[1]))[:10]
print("\nMAIN-SEQUENCE AXIS in module space (top 10 |weight|):")
print("  " + ", ".join(f"{m}({w:+.2f})" for m, w in mods))

# ---- key HR relationships (cosines) ----
def cos(a, b): return float(unit(a) @ unit(b))
L = {lab: M[i] for i, lab in enumerate(labels)}
print("\nKEY RELATIONSHIPS (cosine in 9-PC):")
pairs = [
    ("TRANSFORMATION_EPITHELIAL", "DIFFERENTIATION_COMPLETION"),
    ("TUMOR_TME", "DIFFERENTIATION_COMPLETION"),
    ("TRANSFORMATION_EPITHELIAL", "ACTIVATION_ACUTE"),
    ("NEURODEGEN", "DIFFERENTIATION_COMPLETION"),
    ("AGING_INFLAMMAGING", "ACTIVATION_ACUTE"),
    ("AGING_COLLAPSE", "TRANSFORMATION_EPITHELIAL"),
    ("PULMONARY_FIBROSIS", "LATE_CONSOLIDATION"),
]
for a, b in pairs:
    if a in L and b in L:
        print(f"  cos({a:28s}, {b:26s}) = {cos(L[a], L[b]):+.2f}")
