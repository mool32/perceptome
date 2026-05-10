"""Generate the README hero figure: PC1 × PC4 HPA cell-type scatter with attractor cluster marked.

PC1 = perception breadth (paper 3 Figure 1B axis: myeloid+ to germ−)
PC4 = cancer convergence axis (paper 4.2 Figure 2 axis)

The 8-cell attractor cluster (paper 4.2 P3) is highlighted with larger
markers + black outline + labels — visually demonstrates the central
finding: cancers converge on a non-origin set of "active state"
normal cell types.

Output:
  examples/figures/eigenspace_pc1_pc4.png  (300 dpi, ~50 KB)
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import perceptome as pct  # noqa: E402

OUT_DIR = Path(__file__).resolve().parents[1] / "examples" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "eigenspace_pc1_pc4.png"

# Colorblind-safe palette — 12 distinguishable colors
PALETTE = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
    "#AA3377", "#BBBBBB", "#882255", "#999933", "#CC6677",
    "#117733", "#DDDDDD",
]


def main():
    print("Loading HPA reference + projecting to v0.3 eigenspace...")
    ref = pct.load_hpa_perceptivity()
    R = ref["R"]
    klass = ref["cell_type_class"]

    proj = pct.project(R)
    coords = proj["coordinates"]
    print(f"  projected {len(coords)} cell types into {coords.shape[1]} PCs")

    print("Loading attractor cluster...")
    attr = pct.load_attractor_direction()
    cluster_cells = attr["attractor_cluster_cells"]
    print(f"  attractor cells: {cluster_cells}")
    cluster_present = [c for c in cluster_cells if c in coords.index]
    print(f"  present in HPA: {len(cluster_present)}")

    fig, ax = plt.subplots(figsize=(9, 6.5), dpi=150)

    classes = sorted(klass.unique())
    color_map = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(classes)}

    is_attractor = coords.index.isin(cluster_present)
    for c in classes:
        mask = (klass == c) & ~is_attractor
        ct_in_class = coords.index[mask & (klass == c)]
        if len(ct_in_class) == 0:
            continue
        sub = coords.loc[ct_in_class]
        ax.scatter(
            sub["PC1"], sub["PC4"],
            s=42, alpha=0.65,
            color=color_map[c], edgecolor="white", linewidth=0.5,
            label=c, zorder=2,
        )

    if cluster_present:
        sub_attr = coords.loc[cluster_present]
        ax.scatter(
            sub_attr["PC1"], sub_attr["PC4"],
            s=200, alpha=0.95,
            color="#222222", edgecolor="#FFD700", linewidth=2.0,
            label="cancer attractor cluster (8 cells)", zorder=4, marker="*",
        )

        # Compute cluster centroid + draw a circle to highlight grouping
        cx, cy = sub_attr["PC1"].mean(), sub_attr["PC4"].mean()
        from matplotlib.patches import Circle
        radius = max(
            np.sqrt((sub_attr["PC1"] - cx)**2 + (sub_attr["PC4"] - cy)**2).max() + 0.6,
            1.5,
        )
        circle = Circle((cx, cy), radius, fill=False, ec="#FFD700", lw=1.6, ls="--", alpha=0.7, zorder=3)
        ax.add_patch(circle)

    ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=1)
    ax.axvline(0, color="#bbbbbb", lw=0.5, zorder=1)
    ax.set_xlabel("PC1 — perception breadth\n(myeloid + active states → germ / quiescent)", fontsize=11)
    ax.set_ylabel("PC4 — cancer convergence axis\n(direction tumors take during transformation)", fontsize=11)
    ax.set_title(
        "perceptome v0.2 eigenspace — 154 normal HPA cell types\n"
        "★ = 8-cell attractor cluster that 11 cancers converge toward (Paper 4.2)",
        fontsize=11.5, pad=12,
    )

    # Side legend (cell type classes)
    ax.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5),
        fontsize=8.5, frameon=False,
        title="HPA cell-type class",
    )

    ax.grid(True, alpha=0.18, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"\n  wrote {OUT} ({OUT.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
