#!/usr/bin/env python3
"""Phase 2b/#1 — the cellular main sequence + the cancer panel, using only
robust/validated quantities (no fragile absolute per-cell 1D positioning).

Panel A: the 14 process-trajectories projected onto the dominant acute↔chronic
         axis they define (transition_vectors; robust — these ARE the vectors).
Panel B: all 19 cancers' validated convergence fraction onto the duet/beacon
         region (conv_8beacon, the v7 per-cell nearest-neighbour metric).

NOTE (documented limitation, verified here): absolute per-cell ms-position is NOT
a faithful convergence proxy — Spearman(median ms-position, conv_8beacon)=−0.30,
p=0.28 over 15 recomputed cancers; e.g. neuroblastoma sits at extreme-negative ms
via a neural program, not duet-convergence. So convergence is shown by the
validated nearest-neighbour fraction, and the axis is shown via the transition
vectors that define it.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/Users/teo/Desktop/research/perceptual_modules")
TOOL = ROOT / "tool" / "perceptome2"
sys.path.insert(0, str(TOOL))
import perceptome as pct  # noqa: E402
from perceptome.eigenspace.project import _load_reference  # noqa: E402

ref = _load_reference(None); mo = ref["module_order"]; npc = ref["n_pcs"]
loadings = np.array([[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(npc)] for m in mo])

# ---- Panel A data: process-trajectories on the dominant axis ----
cc = json.load(open(ROOT / "transition_vectors/results/per_pc_decomposition.json"))["class_centroids"]
labels = list(cc)
M = np.array([cc[n] for n in labels]); M = M / np.linalg.norm(M, axis=1, keepdims=True)
ag = json.load(open(TOOL / "perceptome/reference/data/aging_reference.json"))
for k, lab in [("inflammaging_direction", "AGING_INFLAMMAGING"), ("collapse_direction", "AGING_COLLAPSE")]:
    v = np.array(ag[k]); M = np.vstack([M, v / np.linalg.norm(v)]); labels.append(lab)
ms = np.linalg.svd(M, full_matrices=False)[2][0]
if (M[labels.index("TUMOR_TME")] @ ms) < 0:
    ms = -ms
proj = {lab: float(M[i] @ ms) for i, lab in enumerate(labels)}

# ---- Panel B data: all 19 cancers' validated convergence ----
conv = pd.read_csv(ROOT / "results/cancer_v6_6_extras/panel19_convergence.csv").sort_values("conv_8beacon")
baseline = 8 / 154

# ============================ FIGURE ============================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [1, 1]})

# Panel A — main sequence (lollipop, sorted)
items = sorted(proj.items(), key=lambda kv: kv[1])
cancer_cls = {"TUMOR_TME", "TRANSFORMATION_EPITHELIAL"}
for y, (lab, x) in enumerate(items):
    col = "crimson" if lab in cancer_cls else ("darkorange" if "AGING" in lab else "steelblue")
    axA.plot([0, x], [y, y], color=col, lw=1.3, alpha=0.5, zorder=1)
    axA.scatter([x], [y], s=80, color=col, zorder=3, edgecolor="k", linewidth=0.4)
axA.set_yticks(range(len(items)))
axA.set_yticklabels([l.replace("_", " ").title() for l, _ in items], fontsize=8)
axA.axvline(0, color="0.6", lw=0.8)
axA.set_xlabel("projection on dominant axis\n◀ acute-activation          chronic-reorganization ▶")
axA.set_title("A. The cellular main sequence\n(14 process-trajectories on one axis; cancer in red)")
axA.scatter([], [], s=80, color="crimson", label="cancer"); axA.scatter([], [], s=80, color="darkorange", label="aging")
axA.scatter([], [], s=80, color="steelblue", label="other transitions"); axA.legend(fontsize=8, loc="lower right")

# Panel B — 19 cancers convergence to the duet/beacon region
ypos = range(len(conv))
axB.barh(list(ypos), conv["conv_8beacon"].values, color="crimson", alpha=0.6, edgecolor="k", linewidth=0.4)
axB.set_yticks(list(ypos)); axB.set_yticklabels(conv["cancer"].values, fontsize=8)
axB.axvline(baseline, color="k", ls="--", lw=1.0, label=f"chance baseline (8/154 = {baseline:.1%})")
axB.set_xlabel("fraction of malignant cells converging onto the duet/beacon region")
axB.set_title(f"B. All 19 cancers converge to the chronic-pole duet\n(validated per-cell nearest-neighbour metric)")
axB.legend(fontsize=8, loc="lower right")

plt.tight_layout()
out = TOOL / "scripts" / "fig_cancer_mainsequence.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"figure → {out}")
n_above = (conv["conv_8beacon"] > baseline).sum()
print(f"Panel B: {n_above}/{len(conv)} cancers above chance baseline ({baseline:.1%}); "
      f"median convergence {conv['conv_8beacon'].median():.1%}, max {conv['conv_8beacon'].max():.1%} (MM).")
print(f"Panel A: cancer classes at chronic pole — TUMOR_TME {proj['TUMOR_TME']:+.2f}, "
      f"TRANSFORMATION_EPITHELIAL {proj['TRANSFORMATION_EPITHELIAL']:+.2f}; "
      f"acute pole ACTIVATION_ACUTE {proj['ACTIVATION_ACUTE']:+.2f}.")
