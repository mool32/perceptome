"""Generate the README hero figure: PC1 x PC4 HPA cell-type scatter with the
cancer-convergence core marked.

PC1 = perception breadth (paper 3 Figure 1B axis: myeloid+ to germ-)
PC4 = cancer convergence axis

The published convergence result (Spiro 2026, Zenodo 10.5281/zenodo.20542130):
malignant cells of 25 cancers of all lineages converge on the placental
cytotrophoblast -- a normal invasive cell state -- with the megakaryocyte a
lineage-restricted secondary. This figure marks that two-cell core (anchor +
secondary) within the broader "active state" convergence neighbourhood, using
the bundled v0.2 normal-cell eigenspace (the tool's 154 x 44 HPA reference).

Output:
  examples/figures/eigenspace_pc1_pc4.png  (~200 dpi)
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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

ANCHOR = "cytotrophoblasts"        # placental cytotrophoblast — the convergence anchor
SECONDARY = "megakaryocytes"       # lineage-restricted secondary


def main():
    print("Loading HPA reference + projecting to eigenspace...")
    ref = pct.load_hpa_perceptivity()
    R = ref["R"]
    klass = ref["cell_type_class"]
    coords = pct.project(R)["coordinates"]
    print(f"  projected {len(coords)} cell types into {coords.shape[1]} PCs")

    core = [c for c in (ANCHOR, SECONDARY) if c in coords.index]
    if len(core) < 2:
        raise SystemExit(f"core cells not found in HPA reference: have {core}")

    fig, ax = plt.subplots(figsize=(9, 6.5), dpi=150)

    classes = sorted(klass.unique())
    color_map = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(classes)}

    is_core = coords.index.isin(core)
    for c in classes:
        ct = coords.index[(klass == c) & ~is_core]
        if len(ct) == 0:
            continue
        sub = coords.loc[ct]
        ax.scatter(sub["PC1"], sub["PC4"], s=42, alpha=0.6,
                   color=color_map[c], edgecolor="white", linewidth=0.5,
                   label=c, zorder=2)

    # The two-cell core: anchor (large) + secondary (medium)
    a = coords.loc[ANCHOR]
    s = coords.loc[SECONDARY]
    ax.scatter([a["PC1"]], [a["PC4"]], s=430, color="#B4232B", edgecolor="black",
               linewidth=1.6, marker="*", zorder=6)
    ax.scatter([s["PC1"]], [s["PC4"]], s=230, color="#E08A2E", edgecolor="black",
               linewidth=1.3, marker="*", zorder=6)
    ax.annotate("placental cytotrophoblast\n(anchor — 25 cancers converge here)",
                xy=(a["PC1"], a["PC4"]), xytext=(a["PC1"] - 1.0, a["PC4"] + 1.7),
                ha="center", fontsize=9.5, fontweight="bold", color="#B4232B",
                arrowprops=dict(arrowstyle="-", color="#B4232B", lw=1.0))
    ax.annotate("megakaryocyte\n(secondary)",
                xy=(s["PC1"], s["PC4"]), xytext=(s["PC1"] + 1.6, s["PC4"] + 1.6),
                ha="center", fontsize=9, color="#B5651D",
                arrowprops=dict(arrowstyle="-", color="#B5651D", lw=1.0))

    ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=1)
    ax.axvline(0, color="#bbbbbb", lw=0.5, zorder=1)
    ax.set_xlabel("PC1 — perception breadth\n(myeloid + active states -> germ / quiescent)", fontsize=11)
    ax.set_ylabel("PC4 — cancer convergence axis\n(direction tumors take during transformation)", fontsize=11)
    ax.set_title(
        "perceptome eigenspace — 154 normal HPA cell types\n"
        "★ placental cytotrophoblast: the normal-cell anchor that 25 cancers "
        "of all lineages converge on (Spiro 2026)",
        fontsize=11.5, pad=12,
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8.5,
              frameon=False, title="HPA cell-type class")
    ax.grid(True, alpha=0.18, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"  wrote {OUT} ({OUT.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
