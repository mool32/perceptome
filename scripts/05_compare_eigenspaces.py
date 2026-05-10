"""Compare v0.1 (12-PC, 43 modules) vs v0.3 (9-PC, 44 modules) eigenspaces.

Output is the closest-cosine PC mapping table reported in
docs/EIGENSPACE_EVOLUTION.md. Run this script to regenerate the table after
any eigenspace rebuild.
"""

import json
import sys
from pathlib import Path

import numpy as np

V01_REF = Path("/Users/teo/Desktop/research/perceptual_modules/tool/perceptome/perceptome/data/reference_eigenspace.json")
V03_REF = Path(__file__).resolve().parents[1] / "perceptome" / "eigenspace" / "data" / "reference_v03.json"


def loadings_matrix(ref, mods):
    """Build (n_modules × n_pcs) matrix from a reference dict."""
    n_pcs = ref["n_pcs"]
    M = np.zeros((len(mods), n_pcs))
    for j in range(n_pcs):
        loads = ref["loadings"][f"PC{j+1}"]
        for i, m in enumerate(mods):
            M[i, j] = loads.get(m, 0.0)
    return M


def main():
    if not V01_REF.exists():
        print(f"v0.1 reference not found at {V01_REF}; comparison skipped.")
        sys.exit(0)

    with open(V01_REF) as f:
        v01 = json.load(f)
    with open(V03_REF) as f:
        v03 = json.load(f)

    common = sorted(set(v01["module_order"]) & set(v03["module_order"]))
    print(f"Common modules: {len(common)} (v0.1 had {len(v01['module_order'])}, v0.3 has {len(v03['module_order'])})")
    only_v03 = sorted(set(v03["module_order"]) - set(v01["module_order"]))
    print(f"v0.3-only modules: {only_v03}")

    L01 = loadings_matrix(v01, common)
    L03 = loadings_matrix(v03, common)

    # column-normalize
    L01n = L01 / (np.linalg.norm(L01, axis=0, keepdims=True) + 1e-12)
    L03n = L03 / (np.linalg.norm(L03, axis=0, keepdims=True) + 1e-12)

    cos = L01n.T @ L03n  # (v01_pcs × v03_pcs)

    print(f"\nv0.1 PCs: {L01.shape[1]}, v0.3 PCs: {L03.shape[1]}")
    print(f"\nFor each v0.1 PC, closest v0.3 PC (by |cosine| on common modules):\n")
    print(f"{'v0.1':6s} | {'best v0.3':10s} | {'|cos|':>6s} | second-best v0.3 | |cos|")
    print("-" * 72)
    for i in range(L01.shape[1]):
        order = np.argsort(np.abs(cos[i]))[::-1]
        j1, j2 = order[0], order[1]
        print(f"PC{i+1:<4d} | PC{j1+1:<8d} | {abs(cos[i,j1]):.3f} | "
              f"PC{j2+1:<14d} | {abs(cos[i,j2]):.3f}")


if __name__ == "__main__":
    main()
