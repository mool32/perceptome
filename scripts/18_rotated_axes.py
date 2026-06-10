#!/usr/bin/env python3
"""Phase 2b/#2 — define interpretable axes (varimax on PC1-6) and locate the
'main sequence' (chronic-reorganization) axis among them; characterize how the
154 normal cells order along it.
"""
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
mean = np.array(ref["mean_per_module"]); std = np.array(ref["std_per_module"]); std[std < 1e-10] = 1.0
load9 = np.array([[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(npc)] for m in mo])

K = 6  # rotate the bootstrap-stable PC1-6
rot6, Rmat = pct.varimax(load9[:, :K])

R = pct.load_hpa_perceptivity(mode="pseudobulk")["R"].reindex(columns=mo)
cells = list(R.index)
Xz = (R.values.astype(float) - mean) / std
rot_coords = Xz @ rot6                  # 154 × 6 rotated coordinates

print("=" * 76)
print("VARIMAX-ROTATED AXES (PC1-6) — top modules + extreme cells per axis")
print("=" * 76)
for a in range(K):
    w = rot6[:, a]
    top = sorted(zip(mo, w), key=lambda x: -abs(x[1]))[:5]
    c = rot_coords[:, a]; o = np.argsort(c)
    print(f"\nRotated axis R{a+1}:  " + ", ".join(f"{m}({v:+.2f})" for m, v in top))
    print(f"   low:  {[cells[i] for i in o[:4]]}")
    print(f"   high: {[cells[i] for i in o[-4:]]}")

# ---- which rotated axis is the 'main sequence' (chronic-reorganization)? ----
import json
cc = json.load(open(ROOT / "transition_vectors/results/per_pc_decomposition.json"))["class_centroids"]
M = np.array([cc[n] for n in cc]); M = M / np.linalg.norm(M, axis=1, keepdims=True)
ms = np.linalg.svd(M, full_matrices=False)[2][0]   # 9-PC main-sequence axis (#1)
if (M @ ms).sum() < 0:
    ms = -ms

# express each rotated axis in 9-PC: rot6 lives in PC1-6 subspace; pad to 9 via Rmat
# rotated axis a in PC space = Rmat[:,a] (a 6-vector in PC1-6), pad zeros for PC7-9
def cos(a, b): return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
print("\n" + "=" * 76)
print("ALIGNMENT of each rotated axis with the #1 main-sequence (chronic) axis")
print("=" * 76)
for a in range(K):
    rax9 = np.zeros(npc); rax9[:K] = Rmat[:, a]
    print(f"  R{a+1}: cos(rotated-axis, main-sequence) = {cos(rax9, ms):+.2f}")

# ---- best-aligned rotated axis: how do cells order along it? (the static main sequence) ----
aligns = [abs(cos(np.r_[Rmat[:, a], np.zeros(npc - K)], ms)) for a in range(K)]
best = int(np.argmax(aligns))
c = rot_coords[:, best]; o = np.argsort(c)
print(f"\nMAIN-SEQUENCE axis ≈ rotated R{best+1}. Static ordering of normal cells along it:")
print(f"  acute/low pole:    {[cells[i] for i in o[:6]]}")
print(f"  chronic/high pole: {[cells[i] for i in o[-6:]]}")
