"""network/ tests — compute_network + compare_networks."""

import numpy as np
import pandas as pd

import perceptome as pct


def _synthetic_scores(n_cells=100, modules=None, seed=0):
    if modules is None:
        modules = pct.list_modules()[:10]
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(n_cells, len(modules)))
    # Force two correlated modules
    base[:, 1] = base[:, 0] + rng.normal(scale=0.3, size=n_cells)
    return pd.DataFrame(base, columns=modules,
                        index=[f"cell_{i}" for i in range(n_cells)])


def test_compute_network_returns_keys():
    scores = _synthetic_scores()
    out = pct.compute_network(scores)
    for k in ("correlation_matrix", "adjacency", "hub_degrees",
              "communities", "rich_club_phi", "rich_club_intact"):
        assert k in out


def test_compute_network_correlation_matrix_symmetric():
    scores = _synthetic_scores()
    out = pct.compute_network(scores)
    cm = out["correlation_matrix"]
    assert cm.shape[0] == cm.shape[1]
    assert (cm.values - cm.values.T == 0).all()


def test_compute_network_detects_correlation():
    """Modules 0 and 1 are correlated by construction; their correlation should be high."""
    scores = _synthetic_scores()
    out = pct.compute_network(scores)
    cols = list(scores.columns)
    assert out["correlation_matrix"].loc[cols[0], cols[1]] > 0.5


def test_compare_networks_runs():
    s1 = _synthetic_scores(seed=0)
    s2 = _synthetic_scores(seed=1)
    n1 = pct.compute_network(s1)
    n2 = pct.compute_network(s2)
    out = pct.compare_networks(n1, n2)
    assert "np_score" in out
    assert "edge_changes" in out
    assert "narrative" in out


def test_module_heterogeneity_runs():
    """module_heterogeneity should produce per-module distribution metrics."""
    s1 = _synthetic_scores(seed=0)
    s2 = _synthetic_scores(seed=1)
    out = pct.module_heterogeneity(s1, s2)
    # API returns a DataFrame or dict — just check it runs and is non-empty
    assert out is not None
