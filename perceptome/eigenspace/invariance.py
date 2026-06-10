"""Coordinate-invariance check for eigenspace findings.

Turns "let's find better coordinates" (dangerous — invites the tuned-cluster
critique) into "let's show the finding is coordinate-invariant" (bulletproof).
A finding visible in parameter-free PCA AND surviving in an independent nonlinear
geometry (diffusion maps) is geometry, not a method artifact.

Default target: the 8-cell cancer-attractor beacon cluster (Paper 3 / 4.2). Tests
whether the beacon is an unusually TIGHT subset of the 154 normal cell types in
each embedding, vs a null of random equal-size subsets.

All embeddings are fit on the 154 normal HPA cells only. PCA is canonical; varimax
rotation is distance-identical to it (so its compactness equals PCA's by
construction — a built-in confirmation that rotation re-views, not re-geometries);
diffusion maps are the independent nonlinear lens.
"""

import numpy as np

from .project import _load_reference
from .manifold import DiffusionMap, _sqdist
from .rotate import varimax


def _hpa_readiness(module_order, mode="single_cell_scaled"):
    from ..perceptivity.reference import load_hpa_perceptivity
    R = load_hpa_perceptivity(mode=mode)["R"]
    return R.reindex(columns=module_order)


def _zscore(X, mean, std):
    std = np.where(std < 1e-10, 1.0, std)
    return (X - mean) / std


def beacon_compactness(coords, beacon_mask, n_null=10000, seed=0):
    """Mean within-set pairwise distance of the beacon cells vs a null of random
    equal-size subsets. Lower percentile = tighter (more clustered) than chance.

    Returns observed mean pairwise distance, null mean, percentile (∈[0,1]) and z.
    """
    coords = np.asarray(coords, dtype=float)
    n = coords.shape[0]
    idx = np.where(beacon_mask)[0]
    k = len(idx)
    D = np.sqrt(_sqdist(coords))
    tri = np.triu_indices(k, 1)
    obs = D[np.ix_(idx, idx)][tri].mean()

    rng = np.random.default_rng(seed)
    null = np.empty(n_null)
    for b in range(n_null):
        s = rng.choice(n, size=k, replace=False)
        null[b] = D[np.ix_(s, s)][tri].mean()
    return {
        "observed_mean_pairwise": float(obs),
        "null_mean": float(null.mean()),
        "percentile": float((null <= obs).mean()),   # P(random set ≤ beacon) — low = tight
        "z": float((obs - null.mean()) / (null.std() + 1e-12)),
        "k": int(k), "n_null": int(n_null),
    }


def coordinate_invariance(beacon_cells=None, n_components=8, n_null=10000, seed=0,
                          reference=None, mode="single_cell_scaled"):
    """Test whether the beacon cluster is coordinate-invariant across PCA, varimax-
    rotated PCA, and diffusion maps (all fit on the 154 normal HPA cells).

    Parameters
    ----------
    beacon_cells : list[str] | None
        Cell types defining the target cluster. Default = the bundled cancer
        attractor cluster (attractor_v1).
    n_components : int
        Diffusion coordinates to use.
    n_null, seed : int
        Null-distribution controls (shared seed → identical null across embeddings,
        so percentiles are directly comparable).

    Returns
    -------
    dict
        compactness per embedding (pca / rotated_pca / diffusion), `invariant`
        flag, and a `diagnostic` string (PCA-adequate vs curvilinear signal).
    """
    ref = _load_reference(reference)
    module_order = ref["module_order"]
    n_pcs = ref["n_pcs"]
    mean = np.array(ref["mean_per_module"], dtype=float)
    std = np.array(ref["std_per_module"], dtype=float)
    loadings = np.array(
        [[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(n_pcs)] for m in module_order]
    )

    R = _hpa_readiness(module_order, mode=mode)
    cells = list(R.index)
    Z = _zscore(R.values.astype(float), mean, std)

    if beacon_cells is None:
        from ..reference.attractor import load_attractor_direction
        beacon_cells = load_attractor_direction()["attractor_cluster_cells"]
    beacon_set = set(beacon_cells)
    mask = np.array([c in beacon_set for c in cells])

    pca = Z @ loadings
    rotated = Z @ varimax(loadings)[0]
    dm = DiffusionMap(n_components=min(n_components, len(cells) - 2)).fit(Z)

    res = {
        "calibration_mode": mode,   # 'single_cell_scaled' = the beacon finding's native space
        "n_cells": len(cells),
        "beacon_k": int(mask.sum()),
        "beacon_cells_found": sorted(c for c in cells if c in beacon_set),
        "beacon_cells_missing": sorted(beacon_set - set(cells)),
        "diffusion_epsilon": float(dm.epsilon_),
        "space_dims": {"pca": n_pcs, "rotated_pca": n_pcs, "diffusion": dm.n_components},
        "pca": beacon_compactness(pca, mask, n_null, seed),
        "rotated_pca": beacon_compactness(rotated, mask, n_null, seed),
        "diffusion": beacon_compactness(dm.transform(), mask, n_null, seed),
    }

    p_pca = res["pca"]["percentile"]
    p_diff = res["diffusion"]["percentile"]
    res["invariant"] = bool(p_pca < 0.05 and p_diff < 0.05)
    if p_diff < p_pca - 0.02:
        res["diagnostic"] = ("beacon TIGHTER in diffusion than PCA — curvilinear "
                             "manifold signal; worth revisiting in the dynamics epoch")
    elif p_diff <= p_pca + 0.02:
        res["diagnostic"] = ("PCA adequate — diffusion no tighter; the convergence "
                             "is captured by linear geometry")
    else:
        res["diagnostic"] = "beacon LOOSER in diffusion — linear structure dominates"
    return res
