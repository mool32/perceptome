"""Which numbers of the paper does a reference calibration reproduce?

Prints a table of statements from the paper that can be recomputed from a 154 x 44 reference
matrix and the bundled Supporting-Information tables, with the paper's value, the computed value
and a verdict. Run it on the bundled pseudobulk calibration to see where that reference and the
printed numbers part ways, and on the single-cell-scaled calibration (bundled, built from the HPA
table by scripts/13) to see which numbers belong to it:

    python scripts/15_paper_reproduction_report.py                          # bundled pseudobulk
    python scripts/15_paper_reproduction_report.py --calibration single_cell_scaled

Rows come in three groups, following the paper's Methods ("The eigenspace is built on the
pseudobulk HPA reference. For per-cell tumor projection, the reference is re-derived under the
single_cell_scaled calibration"):

  eigenspace (pseudobulk construction)   the PCA of the bundled pseudobulk reference matrix, which is
                                         what fixes the nine components and their loadings
  reference matrix of the calibration    correlations among the 44 modules and the core z-profiles
                                         (Fig. S8, Table S7) of the calibration's own 154 x 44 matrix
  geometry (calibration centroids)       the 154 cell types projected through the fixed loadings under
                                         the calibration: variance shares of the axes (Fig. 1B), the
                                         origins of the 25 cancers, the convergent neighbourhood and the
                                         medoid of the normal manifold (cosine geometry)
"""

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import perceptome as pct  # noqa: E402
from perceptome import convergence as cv  # noqa: E402

PAPER_PAIRS = {frozenset(p) for p in [("UPR-ATF6", "UPR-IRE1"), ("JAK-STAT", "Type I IFN"), ("ERK/MAPK", "NF-κB"),
                                      ("UPR-PERK", "p53"), ("NF-κB", "UPR-PERK"), ("BMP", "Hippo")]}


def _unit(X):
    return X / np.linalg.norm(X, axis=1, keepdims=True)


def rows_for(cal):
    out = []
    add = lambda group, what, paper, got, ok: out.append((group, what, paper, got, "match" if ok else "differs"))

    # -- the eigenspace is built on the pseudobulk reference, whatever the calibration -----------------
    pseudo = pct.load_calibration("pseudobulk")
    pca = pct.rebuild(pseudo.reference_scores[pseudo.module_order].astype(float), bootstrap_n=0)
    ev = np.array(pca["eigenvalues"])
    var = ev / ev.sum()
    g = "eigenspace (pseudobulk construction)"
    add(g, "components with eigenvalue > 1 (Kaiser)", "9", str(int((ev > 1).sum())), (ev > 1).sum() == 9)
    add(g, "variance of the nine components", "73.8%", f"{100 * var[:9].sum():.1f}%", abs(100 * var[:9].sum() - 73.8) < 0.15)

    # -- the calibration's own reference matrix ------------------------------------------------------------
    ref = cal.reference_scores[cal.module_order].astype(float)
    corr = np.corrcoef(ref.values.T)
    iu = np.triu_indices(len(corr), 1)
    big = {frozenset((ref.columns[i], ref.columns[j])) for i, j in zip(*iu) if abs(corr[i, j]) > 0.7}
    s7 = cv.check_calibration_against_table_s7(cal)
    g = "reference matrix of the calibration"
    add(g, "module pairs with |r| > 0.7 of 946 (Fig. S8)", "6", str(len(big)), len(big) == 6)
    add(g, "mean |r| over the 946 module pairs", "0.21", f"{np.abs(corr[iu]).mean():.2f}", abs(np.abs(corr[iu]).mean() - 0.21) < 0.01)
    add(g, "the six pairs are the paper's", "ATF6-IRE1, JAK-STAT-IFN, ERK-NF-kB, PERK-p53, NF-kB-PERK, BMP-Hippo",
        f"{len(big & PAPER_PAIRS)} of 6 shared", big == PAPER_PAIRS)
    add(g, f"Table S7 core z-scores, {s7['n_checked']} printed values (all modules but cGAS-STING)", "equal once rounded",
        f"{s7['n_exact']} of {s7['n_checked']} equal; max |diff| {s7['max_abs_diff_checked']:.4f}", s7["matches"])
    for module, (dm, dc) in s7["known_deviations"].items():
        s7t = cv.load_convergence_reference()["tables"]["S7"].loc[module]
        add(g, f"Table S7, {module} row (flagged in Note S-Q)", f"{s7t['megakaryocyte_z']:+.2f} / {s7t['cytotrophoblast_z']:+.2f}",
            f"{s7t['megakaryocyte_z'] + dm:+.2f} / {s7t['cytotrophoblast_z'] + dc:+.2f}", abs(dm) < 0.006 and abs(dc) < 0.006)

    # -- the 154 cell types in the fixed space under the calibration ------------------------------------
    cent = cal.centroids
    names = list(cent.index)
    S = _unit(cent.values) @ _unit(cent.values).T
    g = "geometry (calibration centroids)"
    v = cent.var(axis=0, ddof=0)
    share = 100 * v / v.sum()
    add(g, "PC1 / PC2 share of variance among the nine axes (Fig. 1B)", "31% / 21%", f"{share.iloc[0]:.1f}% / {share.iloc[1]:.1f}%",
        abs(share.iloc[0] - 31) < 1 and abs(share.iloc[1] - 21) < 1)
    origins = list(cv.load_convergence_reference()["tables"]["S8"]["cell_of_origin"])
    d = [1 - S[names.index(a), names.index(b)] for a, b in itertools.combinations(origins, 2)]
    add(g, "mean pairwise cosine distance among the 25 origins", "0.99", f"{np.mean(d):.3f}", abs(np.mean(d) - 0.99) < 0.006)
    nb = cv.region_members("neighbourhood")
    k = [S[names.index(a), names.index(b)] for a, b in itertools.combinations(nb, 2)]
    rng = np.random.default_rng(0)
    rnd = np.mean([S[i, j] for i, j in (rng.choice(len(names), 2, replace=False) for _ in range(20000))])
    add(g, "four-cell neighbourhood: mutual cosine", "0.77", f"{np.mean(k):.2f}", abs(np.mean(k) - 0.77) < 0.02)
    add(g, "random pairs: mean cosine", "-0.01", f"{rnd:.2f}", abs(rnd + 0.01) < 0.02)
    medoid = names[int((1 - S).sum(axis=1).argmin())]
    add(g, "medoid of the normal manifold (cosine distance, 9-D)", "alveolar cells type 1", medoid, medoid == "alveolar cells type 1")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--calibration", default="pseudobulk", help="'pseudobulk', 'single_cell_scaled' or a calibration JSON")
    args = ap.parse_args(argv)
    cal = pct.load_calibration(args.calibration)
    rows = rows_for(cal)
    df = pd.DataFrame(rows, columns=["group", "statement", "paper", "computed", "verdict"])
    print(f"Calibration: {cal.name}  ({len(cal.reference_scores)} reference cell types)\n")
    print(df.to_markdown(index=False) if _has_tabulate() else df.to_string(index=False))
    n_ok = int((df["verdict"] == "match").sum())
    print(f"\n{n_ok} of {len(df)} rows match the paper.")
    return 0


def _has_tabulate():
    try:
        import tabulate  # noqa: F401
        return True
    except ImportError:
        return False


if __name__ == "__main__":
    sys.exit(main())
