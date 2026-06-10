"""Phase 2 — coordinate-invariance geometry: varimax, diffusion maps, beacon invariance."""

import numpy as np
import pytest

import perceptome as pct


def test_varimax_orthogonal_and_distance_preserving():
    """Varimax is an orthogonal rotation: R is orthonormal and pairwise distances
    in the full PC space are preserved exactly."""
    rng = np.random.default_rng(0)
    loadings = rng.standard_normal((44, 9))
    rotated, R = pct.varimax(loadings)
    assert np.allclose(R @ R.T, np.eye(9), atol=1e-8)        # orthogonal
    X = rng.standard_normal((25, 44))
    d_pca = X @ loadings
    d_rot = X @ rotated
    from scipy.spatial.distance import pdist
    assert np.allclose(pdist(d_pca), pdist(d_rot), atol=1e-8)  # distances preserved


def test_diffusion_oos_round_trip():
    """Nyström out-of-sample reproduces fit coordinates for the fit points."""
    rng = np.random.default_rng(1)
    X = rng.standard_normal((60, 10))
    dm = pct.DiffusionMap(n_components=5).fit(X)
    assert np.allclose(dm.transform_oos(X), dm.transform(), atol=1e-9)


def test_coordinate_invariance_beacon_invariant_scscaled():
    """KEY (Paper 3/4.2): in its native single_cell_scaled calibration the beacon
    cluster is coordinate-invariant — tight in PCA AND diffusion."""
    r = pct.coordinate_invariance(mode="single_cell_scaled", n_null=2000, seed=0)
    assert r["beacon_k"] == 8
    assert r["pca"]["percentile"] < 0.01
    assert r["diffusion"]["percentile"] < 0.01
    assert r["invariant"] is True
    # orthogonal rotation is distance-identical to PCA
    assert r["pca"]["percentile"] == r["rotated_pca"]["percentile"]


def test_coordinate_invariance_positive_control_lineage():
    """Metric validity: a real lineage (immune) is a tight cluster."""
    immune = ["t-cells", "b-cells", "nk-cells", "macrophages",
              "monocytes", "neutrophils", "mast cells", "plasma cells"]
    r = pct.coordinate_invariance(beacon_cells=immune, mode="pseudobulk", n_null=2000, seed=0)
    assert r["pca"]["percentile"] < 0.05


def test_coordinate_invariance_beacon_diffuse_in_pseudobulk():
    """Honest caveat: beacon tightness is calibration-specific — diffuse in the
    pseudobulk (cell-type-level) calibration."""
    r = pct.coordinate_invariance(mode="pseudobulk", n_null=2000, seed=0)
    assert r["pca"]["percentile"] > 0.10          # not tight here
    assert r["invariant"] is False
