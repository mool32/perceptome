"""Generate the README hero figure: the fixed normal-cell space and where 25 cancers land in it.

Panel A  The 154 normal cell types of the bundled reference on the paper's objective coordinate
         plane (Fig. S6A): x = core proximity = (cos -> megakaryocyte + cos -> cytotrophoblast) / 2,
         y = member balance = cos -> megakaryocyte - cos -> cytotrophoblast, cosines taken in the
         standardized nine-dimensional space. Colour = the share of malignant cells whose nearest
         normal cell type this is, as published (Table S2; log scale; grey = below 0.05%).
Panel B  The twelve normal cell types that receive most of the convergence (same Table S2 values),
         against the 1/152 chance baseline; the two-cell core in the accent colour.

Positions come from the bundled single-cell-scaled calibration, the reference the paper projects the
cancers against; colours and bars are the paper's published Table S2 values (computed on the cancer
cells), so positions and colours are not one computation. No cancer data enter the figure.

Output (examples/figures/):
  convergence_region.png, convergence_region_dark.png   (README serves one per colour scheme)
  convergence_region_data.csv                           (the plotted values — the table view)
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import perceptome as pct  # noqa: E402

OUT_DIR = Path(__file__).resolve().parents[1] / "examples" / "figures"

ANCHOR, SECONDARY = "cytotrophoblasts", "megakaryocytes"
CHANCE_PCT = 100 / 152
FLOOR_PCT = 0.05          # below this a cell type is drawn as grey context

# Blue sequential ramp, steps 100..700 of the reference palette
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

# On the dark surface the ramp flips its anchor (more = lighter) and is stepped for that surface
DARK_RAMP = ["#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5", "#5598e7",
             "#6da7ec", "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"]

THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9",
                  axis="#c3c2b7", context="#c3c2b7", ramp=RAMP, accent="#2a78d6", suffix=""),
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781", grid="#2c2c2a",
                 axis="#383835", context="#4d4c48", ramp=DARK_RAMP, accent="#3987e5", suffix="_dark"),
}

LABELS = {
    "cytotrophoblasts": "cytotrophoblast", "megakaryocytes": "megakaryocyte",
    "gastric progenitor cells": "gastric progenitor", "late spermatids": "late spermatid",
    "gastric chief cells": "gastric chief", "pancreatic acinar cells": "pancreatic acinar",
    "esophageal apical cells": "esophageal apical", "hofbauer cells": "Hofbauer cell",
    "migrating cytotrophoblasts": "migrating cytotroph.",
    "epididymal clear cells": "epididymal clear", "differentiating spermatogonia": "diff. spermatogonia",
    "early primary spermatocytes": "early spermatocyte",
}


def plane_table():
    """The 154 normal cell types on the core-proximity plane, with their published occupancy."""
    cent = pct.load_calibration("single_cell_scaled").centroids
    unit = cent.div(np.linalg.norm(cent.values, axis=1), axis=0)
    cos_ct, cos_mk = unit @ unit.loc[ANCHOR], unit @ unit.loc[SECONDARY]
    s2 = pct.load_convergence_reference()["tables"]["S2"]
    return pd.DataFrame({
        "core_proximity": (cos_ct + cos_mk) / 2,
        "member_balance": cos_mk - cos_ct,
        "occupancy_pct": s2["occupancy_pct"].reindex(cent.index),
        "recurrence_of_25": s2["recurrence"].reindex(cent.index),
        "quarantined": s2["qc_status"].reindex(cent.index) == "quarantined",
    }).rename_axis("cell_type")


def draw(df, theme, path):
    t = THEMES[theme]
    cmap = LinearSegmentedColormap.from_list(f"blue_{theme}", t["ramp"])
    norm = LogNorm(vmin=FLOOR_PCT, vmax=14)
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": t["ink"], "axes.labelcolor": t["ink2"],
                         "xtick.color": t["muted"], "ytick.color": t["muted"], "axes.edgecolor": t["axis"]})

    fig = plt.figure(figsize=(11.6, 6.3), dpi=200, facecolor=t["surface"])
    ax = fig.add_axes([0.068, 0.345, 0.50, 0.505])          # panel A: the map
    cax = fig.add_axes([0.125, 0.222, 0.34, 0.020])         # its horizontal colour scale
    bx = fig.add_axes([0.765, 0.345, 0.22, 0.505])          # panel B: the ranking
    for a in (ax, bx):
        a.set_facecolor(t["surface"])
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            a.spines[side].set_color(t["axis"])
        a.grid(True, color=t["grid"], linewidth=0.8, zorder=0)
        a.set_axisbelow(True)

    # -- panel A: the map --
    ring = dict(edgecolors=t["surface"], linewidths=1.1)
    quar = df[df["quarantined"]]
    ctx = df[~df["quarantined"] & ~(df["occupancy_pct"] >= FLOOR_PCT)]
    hot = df[~df["quarantined"] & (df["occupancy_pct"] >= FLOOR_PCT)].sort_values("occupancy_pct")
    ax.scatter(ctx["core_proximity"], ctx["member_balance"], s=46, c=t["context"], zorder=2, **ring)
    ax.scatter(quar["core_proximity"], quar["member_balance"], s=46, facecolors="none",
               edgecolors=t["muted"], linewidths=1.3, zorder=2)
    sc = ax.scatter(hot["core_proximity"], hot["member_balance"], s=60, c=hot["occupancy_pct"], cmap=cmap,
                    norm=norm, zorder=3, **ring)
    core = df.loc[[ANCHOR, SECONDARY]]
    ax.scatter(core["core_proximity"], core["member_balance"], s=330, marker="*", c=core["occupancy_pct"],
               cmap=cmap, norm=norm, zorder=4, edgecolors=t["ink"], linewidths=1.0)

    offsets = {  # label offsets in points (dx, dy); thin leader lines keep each label attached to its mark
        ANCHOR: (-20, -34), SECONDARY: (-26, 40), "gastric progenitor cells": (-60, -30),
        "late spermatids": (-74, 16),
    }
    for name, (dx, dy) in offsets.items():
        r = df.loc[name]
        weight = "bold" if name in (ANCHOR, SECONDARY) else "normal"
        txt = LABELS[name] + {ANCHOR: "  (anchor)", SECONDARY: "  (secondary)"}.get(name, "")
        ax.annotate(txt, (r["core_proximity"], r["member_balance"]), xytext=(dx, dy),
                    textcoords="offset points", fontsize=10.5, color=t["ink"], fontweight=weight,
                    ha="center", va="center", zorder=5,
                    arrowprops=dict(arrowstyle="-", color=t["muted"], linewidth=0.9, shrinkA=0, shrinkB=5))
    ax.set_xlim(-0.88, 0.95)
    ax.set_ylim(-0.86, 0.80)
    ax.set_xlabel("core proximity   (mean cosine to the cytotrophoblast and the megakaryocyte)", fontsize=10.5)
    ax.set_ylabel("member balance   (\u2190 cytotrophoblast  \u00b7  megakaryocyte \u2192)", fontsize=10.5)
    ax.tick_params(labelsize=10)
    ax.set_title("A   154 normal cell types in the fixed reference space", loc="left", fontsize=12.5,
                 fontweight="bold", color=t["ink"], pad=10)

    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_ticks([0.1, 1, 10])
    cb.set_ticklabels(["0.1", "1", "10"])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=10, colors=t["muted"], length=3)
    cb.outline.set_visible(False)
    cb.set_label("% of malignant cells (25 cancers) whose nearest normal cell type this is   (log scale; grey: < 0.05%)",
                 fontsize=9.8, color=t["ink2"], labelpad=6, loc="left")

    # -- panel B: the ranking --
    top = df.dropna(subset=["occupancy_pct"]).sort_values("occupancy_pct", ascending=False).head(12)[::-1]
    colors = [t["accent"] if n in (ANCHOR, SECONDARY) else t["context"] for n in top.index]
    y = np.arange(len(top))
    bx.barh(y, top["occupancy_pct"], height=0.62, color=colors, zorder=3)
    bx.axvline(CHANCE_PCT, color=t["ink2"], linewidth=1.2, zorder=4)
    bx.text(CHANCE_PCT + 0.2, -0.95, "chance 0.66%", fontsize=9.5, color=t["ink2"], va="center", ha="left")
    bx.set_yticks(y)
    bx.set_yticklabels([LABELS[n] for n in top.index], fontsize=10.5, color=t["ink2"])
    bx.tick_params(axis="y", length=0)
    bx.tick_params(axis="x", labelsize=10)
    bx.set_xlim(0, 15.8)
    bx.set_ylim(-1.3, len(top) - 0.4)
    bx.grid(False, axis="y")
    for n in (ANCHOR, SECONDARY):
        i = list(top.index).index(n)
        bx.text(top.loc[n, "occupancy_pct"] + 0.3, i, f"{top.loc[n, 'occupancy_pct']:.1f}%", va="center",
                fontsize=10.5, fontweight="bold", color=t["ink"])
    bx.set_xlabel("% of malignant cells (Table S2)", fontsize=10.5)
    bx.set_title("B   Where the cancers land", loc="left", fontsize=12.5, fontweight="bold",
                 color=t["ink"], pad=10)

    fig.text(0.068, 0.945, "Diverse cancers converge on one region of normal cell-state space",
             fontsize=13.5, fontweight="bold", color=t["ink"], ha="left", va="center")
    fig.text(0.068, 0.070,
             "Positions: bundled single-cell-scaled calibration (standardized 9-D cosines).  Colours and bars: published results "
             "(Spiro 2026, Table S2).\n"
             "Open circles: erythrocyte and platelet, quarantined low-complexity artifacts.\n"
             "Core = cytotrophoblast (anchor) + megakaryocyte, which fails the transcriptome-complexity control.",
             fontsize=9.3, color=t["ink2"], ha="left", va="center", linespacing=1.55)
    fig.savefig(path, dpi=200, facecolor=t["surface"])
    plt.close(fig)
    print(f"  wrote {path} ({path.stat().st_size / 1024:.0f} KB)")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = plane_table()
    for theme, t in THEMES.items():
        draw(df, theme, OUT_DIR / f"convergence_region{t['suffix']}.png")
    csv = OUT_DIR / "convergence_region_data.csv"
    df.round(4).to_csv(csv)
    print(f"  wrote {csv} ({len(df)} rows)")


if __name__ == "__main__":
    main()
