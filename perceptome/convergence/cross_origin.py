"""Cross-origin contraction of malignant states (paper Fig. 2, Methods, Note S-G).

Each cancer's malignant centroid in the nine-dimensional reference space (the mean over its
quality-controlled malignant cells) is paired with its normal cell of origin (Table S8). Unrelated
normal origins are near-orthogonal (mean pairwise cosine distance 0.99; three of the 300 origin
pairs share an origin and contribute exactly 0); the malignant centroids are not (about 0.16).

The paper then takes its own headline apart. The raw contraction is largely a matched-averaging and
shared-magnitude effect: aggregates of unrelated normal cells formed by the same averaging reach
about 0.13. The fold-contraction against random draws of 25 normal cell types (a Stayton-type test)
mis-specified the null and is withdrawn. What survives is a small convergence in module shape, which
the paper reports as a specification curve rather than one number. This module computes the
quantities of that analysis — ``contraction`` (the observed distances and their ratio to the origins')
and ``contraction_null`` (the matched-averaging null in generic form) — and leaves the specification
to the caller: run it on coordinates from ``assign_cells(..., shape_only=True)``, or on module scores,
and report the variants side by side.
"""

import numpy as np
import pandas as pd

from ..eigenspace.calibration import Calibration, load_calibration


def _distance_matrix(A):
    A = np.asarray(A, dtype=float)
    norm = np.linalg.norm(A, axis=1, keepdims=True)
    U = A / np.where(norm < 1e-12, np.nan, norm)
    D = 1.0 - U @ U.T
    np.fill_diagonal(D, 0.0)
    return D


def _mean_pairwise(A):
    D = _distance_matrix(A)
    iu = np.triu_indices(D.shape[0], 1)
    return float(np.nanmean(D[iu])) if len(iu[0]) else float("nan")


def cosine_distances(X):
    """Pairwise cosine distance (1 - cosine) between the rows of ``X``, as a DataFrame."""
    idx = X.index if hasattr(X, "index") else None
    return pd.DataFrame(_distance_matrix(X), index=idx, columns=idx)


def mean_pairwise_cosine_distance(X):
    """Mean cosine distance over the n(n-1)/2 pairs of rows of ``X``.

    0 for identical directions, 1 for orthogonal ones, 2 for opposite ones: unrelated states of a
    high-dimensional space sit near 1.
    """
    return _mean_pairwise(X)


def contraction(coordinates, group, origins=None, calibration="pseudobulk", shape_only=False,
                min_cells=1):
    """Mean pairwise cosine distance among group centroids, and among their normal origins.

    Parameters
    ----------
    coordinates : DataFrame
        Cells x dimensions: the standardized nine-dimensional coordinates returned by
        ``assign_cells(...)['coordinates']`` (or any other cell-level matrix, such as module
        scores, if no ``origins`` are given).
    group : Series
        Cancer / cohort label per cell (aligned to ``coordinates`` by index).
    origins : Series | dict | None
        Normal cell type of origin of each group (names of the reference cell types; Table S8 of
        the paper maps its 25 cancers). If given, the same distance is computed among the origins'
        reference centroids — one row per group, repeated origins repeated, so groups sharing an
        origin contribute pairs at distance 0, as in the paper.
    calibration : 'pseudobulk' | 'single_cell_scaled' | Calibration
        The calibration whose reference centroids stand for the origins; use the one the cells
        were assigned with.
    shape_only : bool
        Take the origins from the shape-only variant of the calibration (use together with
        ``assign_cells(..., shape_only=True)``).
    min_cells : int
        Drop groups with fewer cells than this.

    Returns
    -------
    dict
        centroids, distances      group centroids and their pairwise cosine distance matrix
        mean_distance             mean over the n(n-1)/2 group pairs (paper: ~0.16 for 25 cancers)
        n_groups, n_cells         number of groups, cells per group
        origin_distances, origin_mean_distance, ratio, n_origin_pairs_shared
                                  (with ``origins``) the same for the origins (paper: ~0.99); ratio =
                                  origin_mean_distance / mean_distance; number of pairs of groups
                                  that share an origin
    """
    X = pd.DataFrame(coordinates)
    g = pd.Series(group).reindex(X.index)
    keep = g.notna()
    X, g = X[keep], g[keep]
    counts = g.value_counts()
    kept = counts[counts >= min_cells].index
    X, g = X[g.isin(kept)], g[g.isin(kept)]
    if g.nunique() < 2:
        raise ValueError("contraction needs at least two groups with enough cells.")

    centroids = X.groupby(g, sort=True).mean()
    out = {
        "centroids": centroids,
        "distances": cosine_distances(centroids),
        "mean_distance": _mean_pairwise(centroids.values),
        "n_groups": int(len(centroids)),
        "n_cells": counts.reindex(centroids.index).astype(int),
    }
    if origins is not None:
        o = pd.Series(origins)
        absent = [k for k in centroids.index if k not in o.index]
        if absent:
            raise ValueError(f"origins has no entry for: {absent}")
        cal = calibration if isinstance(calibration, Calibration) else load_calibration(calibration)
        if shape_only:
            cal = cal.shape_only_calibration()
        ref = cal.centroids
        if ref.shape[1] != centroids.shape[1]:
            raise ValueError(f"origins live in the {ref.shape[1]}-dimensional reference space but the "
                             f"coordinates have {centroids.shape[1]} columns; pass assign_cells(...)['coordinates'].")
        names = o.loc[centroids.index]
        unknown = sorted(set(names) - set(ref.index))
        if unknown:
            raise ValueError(f"origins not among the reference cell types: {unknown}")
        oc = ref.loc[names.values].copy()
        oc.index = centroids.index
        out["origin_distances"] = cosine_distances(oc)
        out["origin_mean_distance"] = _mean_pairwise(oc.values)
        out["ratio"] = float(out["origin_mean_distance"] / out["mean_distance"]) if out["mean_distance"] else float("inf")
        n = len(names)
        out["n_origin_pairs_shared"] = int(sum(names.iloc[i] == names.iloc[j]
                                               for i in range(n) for j in range(i + 1, n)))
    return out


def contraction_null(coordinates, population, sizes, observed=None, n_draws=300, seed=0):
    """Matched-averaging null of the contraction (Note S-G), in generic form.

    Averaging many cells shrinks them toward a shared centre, so tight centroids are expected even
    from unrelated states. The null forms aggregates of normal cells by the same averaging and of the
    same sizes, and asks how tight *their* centroids are: per draw, one aggregate per entry of
    ``sizes`` — ``size`` cells of one normal micro-population, averaged — and the mean pairwise cosine
    distance among the aggregates. Distinct micro-populations are used within a draw when there are
    enough of them. The paper's null (200-300 draws over a seven-tissue pool of coherent normal
    micro-populations; mean 0.129, 95% interval [0.074, 0.213] for its 25 cancers) is built from
    normal single-cell data that are not part of this package: pass your own normal cells, projected
    with the same calibration and the same ``shape_only`` setting as the cancers.

    Parameters
    ----------
    coordinates : DataFrame
        Normal cells x dimensions (``assign_cells(...)['coordinates']`` of normal cells).
    population : Series
        Micro-population label per normal cell (for example author cell type within tissue). What
        counts as coherent is the caller's definition.
    sizes : sequence of int
        Aggregate sizes — typically each cancer's number of malignant cells, or one common size.
    observed : float | None
        The observed mean pairwise distance (``contraction(...)['mean_distance']``).
    n_draws, seed : int
        Number of null draws and the seed of the generator.

    Returns
    -------
    dict
        null           array of ``n_draws`` null values
        null_mean, interval95   mean and 2.5-97.5 percentiles of the null
        With ``observed``: n_below (draws below it), percentile (their share, %), p_lower
        = (1 + #{null <= observed}) / (n_draws + 1), the one-sided permutation P the paper uses
        ((8+1)/(300+1) = 0.030 for its shape-only result).
    """
    X = np.asarray(pd.DataFrame(coordinates), dtype=float)
    pop = pd.Series(population, index=pd.DataFrame(coordinates).index)
    members = {p: np.flatnonzero((pop == p).values) for p in pop.dropna().unique()}
    names = sorted(members, key=str)
    sizes = [int(s) for s in sizes]
    if len(names) < 1 or len(sizes) < 2:
        raise ValueError("the null needs at least one micro-population and two aggregate sizes.")
    if min(sizes) < 1:
        raise ValueError("aggregate sizes must be positive.")
    distinct = len(names) >= len(sizes)

    rng = np.random.default_rng(seed)
    null = np.empty(int(n_draws))
    for t in range(len(null)):
        chosen = rng.choice(len(names), size=len(sizes), replace=not distinct)
        cents = np.empty((len(sizes), X.shape[1]))
        for a, (ci, size) in enumerate(zip(chosen, sizes)):
            idx = members[names[ci]]
            take = rng.choice(idx, size=size, replace=len(idx) < size)
            cents[a] = X[take].mean(axis=0)
        null[t] = _mean_pairwise(cents)

    out = {"null": null, "null_mean": float(np.nanmean(null)),
           "interval95": (float(np.nanpercentile(null, 2.5)), float(np.nanpercentile(null, 97.5))),
           "n_draws": int(len(null))}
    if observed is not None:
        out["observed"] = float(observed)
        out["n_below"] = int((null < observed).sum())
        out["percentile"] = float(100 * out["n_below"] / len(null))
        out["p_lower"] = float((1 + (null <= observed).sum()) / (len(null) + 1))
    return out
