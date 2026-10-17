"""Cross-origin contraction (Fig. 2, Methods, Note S-G): distances, origin pairing, matched-averaging null."""

import itertools

import numpy as np
import pandas as pd
import pytest

import perceptome as pct
from perceptome import convergence as cv


@pytest.fixture(scope="module")
def origins():
    return cv.load_convergence_reference()["tables"]["S8"]["cell_of_origin"]


@pytest.fixture(scope="module")
def centroids():
    return pct.load_calibration("pseudobulk").centroids


def cohort(origins, centroids, alpha, n=40, noise=0.3, seed=0):
    """Cells of each cancer: its origin's reference centroid pulled a fraction ``alpha`` toward the core."""
    rng = np.random.default_rng(seed)
    core = centroids.loc[["cytotrophoblasts", "megakaryocytes"]].mean()
    rows, groups = [], []
    for cancer, o in origins.items():
        c = (1 - alpha) * centroids.loc[o] + alpha * core
        rows.append(c.values + rng.normal(0, noise, (n, centroids.shape[1])))
        groups += [cancer] * n
    idx = [f"c{i}" for i in range(len(groups))]
    return pd.DataFrame(np.vstack(rows), index=idx, columns=centroids.columns), pd.Series(groups, index=idx)


def test_cosine_distances_are_zero_one_two_for_parallel_orthogonal_opposite():
    X = pd.DataFrame([[1, 0], [0, 1], [-1, 0], [2, 0]], index=list("abcd"), dtype=float)
    D = cv.cosine_distances(X)
    assert D.loc["a", "b"] == pytest.approx(1.0)
    assert D.loc["a", "c"] == pytest.approx(2.0)
    assert D.loc["a", "d"] == pytest.approx(0.0)
    np.testing.assert_allclose(np.diag(D), 0.0)
    assert cv.mean_pairwise_cosine_distance(X.iloc[:2]) == pytest.approx(1.0)
    assert cv.mean_pairwise_cosine_distance(X.iloc[:3]) == pytest.approx((1 + 2 + 1) / 3)


def test_table_s8_has_300_pairs_of_which_three_share_an_origin(origins):
    """Results: 'three of the 300 origin pairs share an assigned origin by construction and contribute exactly 0'."""
    pairs = list(itertools.combinations(origins.index, 2))
    shared = [(a, b) for a, b in pairs if origins[a] == origins[b]]
    assert len(origins) == 25 and len(pairs) == 300
    assert set(shared) == {("B-ALL", "DLBCL"), ("Germ-cell", "Seminoma"), ("Osteosarcoma", "Ewing sarcoma")}


def test_the_origins_of_the_25_cancers_are_near_orthogonal(origins, centroids):
    """The paper prints 0.99 (single-cell-scaled calibration); the bundled pseudobulk centroids give 0.976."""
    X, g = cohort(origins, centroids, alpha=0.0, noise=0.0, n=2)
    r = cv.contraction(X, g, origins=origins)
    assert r["n_groups"] == 25 and r["n_origin_pairs_shared"] == 3
    assert r["origin_mean_distance"] == pytest.approx(0.976, abs=0.003)
    assert r["mean_distance"] == pytest.approx(r["origin_mean_distance"], abs=1e-9)
    assert r["ratio"] == pytest.approx(1.0)
    D = r["origin_distances"]
    assert D.loc["B-ALL", "DLBCL"] == pytest.approx(0.0, abs=1e-12)
    assert D.loc["HCC", "BRCA"] > 0.5


def test_cancers_pulled_toward_a_common_region_contract(origins, centroids):
    X, g = cohort(origins, centroids, alpha=0.9)
    r = cv.contraction(X, g, origins=origins)
    assert r["mean_distance"] < 0.05
    assert r["ratio"] > 20
    half = cv.contraction(*cohort(origins, centroids, alpha=0.5), origins=origins)
    assert r["mean_distance"] < half["mean_distance"] < half["origin_mean_distance"]


def test_contraction_without_origins_and_with_module_scores():
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(0, 1, (60, 44)), columns=pct.list_modules())
    g = pd.Series(np.repeat(list("abcdef"), 10), index=X.index)
    r = cv.contraction(X, g)
    assert r["n_groups"] == 6 and "ratio" not in r
    assert r["distances"].shape == (6, 6) and list(r["n_cells"]) == [10] * 6


def test_contraction_applies_the_cell_floor_and_validates_inputs(origins, centroids):
    X, g = cohort(origins, centroids, alpha=0.2, n=12)
    g.iloc[:6] = "tiny"
    kept = cv.contraction(X, g, min_cells=10)
    assert "tiny" not in kept["centroids"].index
    with pytest.raises(ValueError, match="at least two groups"):
        cv.contraction(X, pd.Series("one", index=X.index))
    with pytest.raises(ValueError, match="no entry"):
        cv.contraction(X, g, origins={"HCC": "hepatocytes"})
    with pytest.raises(ValueError, match="not among the reference"):
        cv.contraction(X, g, origins={**origins.to_dict(), "tiny": "no such cell"})
    idx = [f"m{i}" for i in range(40)]
    modules = pd.DataFrame(np.random.default_rng(0).normal(size=(40, 44)), index=idx)
    with pytest.raises(ValueError, match="dimensional"):                 # origins live in the 9-D space
        cv.contraction(modules, pd.Series(["HCC", "BRCA"] * 20, index=idx), origins=origins)


def test_shape_only_origins_come_from_the_shape_only_calibration(origins, centroids):
    X, g = cohort(origins, centroids, alpha=0.0, noise=0.0, n=2)
    plain = cv.contraction(X, g, origins=origins)
    shape = cv.contraction(X, g, origins=origins, shape_only=True)
    assert not np.isclose(plain["origin_mean_distance"], shape["origin_mean_distance"])


# ── matched-averaging null ───────────────────────────────────────────────────

def normal_cells(n_pop=12, n_cells=30, dim=9, shared=0.0, noise=0.5, seed=0):
    rng = np.random.default_rng(seed)
    means = rng.normal(0, 1, (n_pop, dim))
    offset = np.zeros(dim)
    offset[0] = shared
    X = np.vstack([m + offset + rng.normal(0, noise, (n_cells, dim)) for m in means])
    idx = [f"n{i}" for i in range(len(X))]
    return pd.DataFrame(X, index=idx), pd.Series(np.repeat(np.arange(n_pop), n_cells), index=idx)


def test_the_null_reports_the_papers_permutation_p_and_is_reproducible():
    X, pop = normal_cells()
    a = cv.contraction_null(X, pop, [20] * 6, observed=0.0, n_draws=300, seed=4)
    b = cv.contraction_null(X, pop, [20] * 6, observed=0.0, n_draws=300, seed=4)
    np.testing.assert_array_equal(a["null"], b["null"])
    assert a["n_draws"] == 300 and a["n_below"] == 0 and a["percentile"] == 0.0
    assert a["p_lower"] == pytest.approx(1 / 301)                  # (0 + 1) / (300 + 1)
    assert a["interval95"][0] < a["null_mean"] < a["interval95"][1]
    hi = cv.contraction_null(X, pop, [20] * 6, observed=np.inf, n_draws=300, seed=4)
    assert hi["n_below"] == 300 and hi["p_lower"] == pytest.approx(1.0)
    other = cv.contraction_null(X, pop, [20] * 6, n_draws=300, seed=5)
    assert not np.array_equal(a["null"], other["null"]) and "p_lower" not in other


def test_the_papers_permutation_p_arithmetic():
    """Note S-G (ii): 8 of 300 draws at or below the observed value, P = (8+1)/(300+1) = 0.030."""
    null = np.concatenate([np.full(8, 0.01), np.full(292, 0.05)])
    observed = 0.02
    assert (1 + (null <= observed).sum()) / (len(null) + 1) == pytest.approx(0.030, abs=0.0005)
    assert (null > observed).sum() == 292


def test_a_shared_magnitude_makes_unrelated_states_look_contracted():
    """The paper's point: contraction of unrelated states is expected when they share a large component."""
    loose, pop = normal_cells(shared=0.0)
    tight, _ = normal_cells(shared=12.0)
    a = cv.contraction_null(loose, pop, [20] * 6, n_draws=100, seed=1)["null_mean"]
    b = cv.contraction_null(tight, pop, [20] * 6, n_draws=100, seed=1)["null_mean"]
    assert a > 0.5 > 0.1 > b


def test_the_null_reuses_populations_when_there_are_too_few_and_validates():
    X, pop = normal_cells(n_pop=3, n_cells=8)
    out = cv.contraction_null(X, pop, [20] * 6, n_draws=20, seed=0)      # 6 aggregates, 3 populations, 8 cells each
    assert len(out["null"]) == 20 and np.isfinite(out["null"]).all()
    with pytest.raises(ValueError):
        cv.contraction_null(X, pop, [20], n_draws=5)
    with pytest.raises(ValueError):
        cv.contraction_null(X, pop, [20, 0], n_draws=5)
