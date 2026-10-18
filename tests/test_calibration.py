"""Standardized projection + reference calibrations (the paper's coordinate system)."""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import perceptome as pct
from perceptome.eigenspace.project import _load_reference


@pytest.fixture(scope="module")
def ref():
    return _load_reference()


@pytest.fixture(scope="module")
def hpa_R():
    # the bundled matrix is float32; do the arithmetic in float64 so tolerances stay tight
    return pct.load_hpa_perceptivity()["R"].astype(float)


def _loadings(ref):
    mods = ref["module_order"]
    return np.array([[ref["loadings"][f"PC{j + 1}"][m] for j in range(ref["n_pcs"])] for m in mods])


# ── project(): z-score against the reference, then the fixed loadings ─────────

def test_project_is_standardized_by_default(ref, hpa_R):
    mods = ref["module_order"]
    z = (hpa_R[mods].values - np.array(ref["mean_per_module"])) / np.array(ref["std_per_module"])
    expected = z @ _loadings(ref)
    out = pct.project(hpa_R)
    np.testing.assert_allclose(out["coordinates"].values, expected, atol=1e-8)
    assert out["calibration"] == "pseudobulk"


def test_reference_cell_types_are_centered_with_eigenvalue_variance(ref, hpa_R):
    coords = pct.project(hpa_R)["coordinates"]
    np.testing.assert_allclose(coords.mean(axis=0).values, 0.0, atol=1e-6)
    np.testing.assert_allclose(coords.var(axis=0, ddof=0).values, ref["eigenvalues"][:9], rtol=1e-6)


def test_legacy_projection_is_still_available(ref, hpa_R):
    mods = ref["module_order"]
    out = pct.project(hpa_R, standardize=False)
    np.testing.assert_allclose(out["coordinates"].values, hpa_R[mods].values @ _loadings(ref), atol=1e-6)
    assert out["calibration"] is None


def _mean_pairwise_cosine(X):
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
    iu = np.triu_indices(len(X), 1)
    return float((Xn @ Xn.T)[iu].mean())


def test_normal_cell_types_are_near_orthogonal_only_when_standardized(hpa_R):
    """Paper: states of unrelated origin are near-orthogonal (mean cosine distance ~0.99).

    The unscaled projection shared one large offset (mean cosine ~0.97) — the defect fixed here.
    """
    std = pct.project(hpa_R)["coordinates"].values
    legacy = pct.project(hpa_R, standardize=False)["coordinates"].values
    assert abs(_mean_pairwise_cosine(std)) < 0.02
    assert _mean_pairwise_cosine(legacy) > 0.9


def test_pc1_runs_from_inflammatory_to_terminal_germline(hpa_R):
    """Paper Results: PC1 'runs from inflammatory cells to terminal germline cells'."""
    pc1 = pct.project(hpa_R)["coordinates"]["PC1"]
    inflammatory = ["neutrophils", "monocytes", "macrophages"]
    germline = ["oocytes", "megakaryocytes"]
    assert (pc1[inflammatory] < 0).all()
    assert (pc1[germline] > 0).all()
    assert pc1[inflammatory].max() < pc1[germline].min()


def test_missing_modules_sit_at_the_reference_mean(ref, hpa_R):
    sample = hpa_R.iloc[:5]
    dropped = pct.project(sample.drop(columns=["NPAS4"]))["coordinates"]
    neutral = sample.copy()
    neutral["NPAS4"] = ref["mean_per_module"][ref["module_order"].index("NPAS4")]
    np.testing.assert_allclose(dropped.values, pct.project(neutral)["coordinates"].values, atol=1e-8)


def test_missing_modules_argument_neutralizes_uncovered_modules(hpa_R):
    sample = hpa_R.iloc[:3].copy()
    zeroed = sample.copy()
    zeroed["NPAS4"] = 0.0   # what score_modules writes for a module without genes
    a = pct.project(zeroed, missing_modules=["NPAS4"], warn_scale=False)["coordinates"]
    b = pct.project(sample.drop(columns=["NPAS4"]), warn_scale=False)["coordinates"]
    np.testing.assert_allclose(a.values, b.values, atol=1e-8)


# ── shifts ────────────────────────────────────────────────────────────────────

def test_project_shift_equals_difference_of_projections(hpa_R):
    a, b = hpa_R.iloc[:4], hpa_R.iloc[10:14]
    delta = pd.DataFrame(a.values - b.values, columns=hpa_R.columns)
    shift = pct.project_shift(delta)
    diff = pct.project(a)["coordinates"].values - pct.project(b)["coordinates"].values
    np.testing.assert_allclose(shift.values, diff, atol=1e-8)


def test_project_shift_accepts_a_series(hpa_R):
    delta = hpa_R.iloc[0] - hpa_R.iloc[1]
    v = pct.project_shift(delta)
    assert isinstance(v, pd.Series) and len(v) == 9


def test_compare_conditions_delta_eigenspace_is_the_difference_of_projected_means():
    rng = np.random.default_rng(1)
    mods = pct.list_modules()
    s1 = pd.DataFrame(rng.normal(2.5, 0.4, size=(80, 44)), columns=mods)
    s2 = s1 + rng.normal(0.3, 0.1, size=44)
    out = pct.compare_conditions(s2, s1, project_to_eigenspace=True)
    m2 = pct.project(s2.mean().to_frame().T, warn_scale=False)["coordinates"].iloc[0].values
    m1 = pct.project(s1.mean().to_frame().T, warn_scale=False)["coordinates"].iloc[0].values
    np.testing.assert_allclose(out["delta_eigenspace"].values, m2 - m1, atol=1e-8)


# ── scale diagnostic ─────────────────────────────────────────────────────────

def test_per_cell_scale_scores_trigger_the_calibration_warning(hpa_R):
    with pytest.warns(UserWarning, match="single_cell_scaled"):
        pct.project(hpa_R.iloc[:20] / 10.0)


def test_reference_scale_scores_do_not_warn(hpa_R):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        pct.project(hpa_R.iloc[:20])


# ── Calibration object ───────────────────────────────────────────────────────

def test_bundled_calibration_matches_reference_statistics(ref):
    cal = pct.load_calibration("pseudobulk")
    np.testing.assert_allclose(cal.mean.values, ref["mean_per_module"], atol=1e-8)
    np.testing.assert_allclose(cal.std.values, ref["std_per_module"], atol=1e-8)
    assert cal.centroids.shape == (154, 9)
    assert "cytotrophoblasts" in cal.centroids.index


def test_calibration_dict_roundtrip(tmp_path):
    cal = pct.load_calibration("pseudobulk")
    path = tmp_path / "cal.json"
    cal.save(path)
    back = pct.load_calibration(str(path))
    np.testing.assert_allclose(back.centroids.values, cal.centroids.values, atol=1e-10)
    assert back.name == "pseudobulk"


def test_shape_only_is_invariant_to_a_constant_shift_of_every_module(hpa_R):
    cal = pct.load_calibration("pseudobulk").shape_only_calibration()
    sample = hpa_R.iloc[:6]
    a = cal.project(sample).values
    b = cal.project(sample + 1.7).values
    np.testing.assert_allclose(a, b, atol=1e-8)
    plain = pct.load_calibration("pseudobulk")
    assert not np.allclose(plain.project(sample).values, plain.project(sample + 1.7).values)


def test_shape_only_reference_is_centred_per_cell_type():
    cal = pct.load_calibration("pseudobulk").shape_only_calibration()
    np.testing.assert_allclose(cal.reference_scores.mean(axis=1).values, 0.0, atol=1e-8)
    assert cal.name.endswith("+shape_only")


def test_single_cell_scaled_calibration_must_be_built_from_hpa():
    path = Path(pct.__file__).parent / "eigenspace" / "data" / "calibration_single_cell_scaled_v03.json"
    if path.exists():
        cal = pct.load_calibration("single_cell_scaled")
        assert cal.reference_scores.shape == (154, 44)
    else:
        with pytest.raises(FileNotFoundError, match="13_build_single_cell_calibration"):
            pct.load_calibration("single_cell_scaled")


def test_unknown_calibration_path_raises():
    with pytest.raises(FileNotFoundError):
        pct.load_calibration("/nonexistent/cal.json")


# ── hpa_reference_scores ─────────────────────────────────────────────────────

def _toy_expression(n_genes=6000):
    rng = np.random.default_rng(0)
    genes = pct.get_genes("NF-κB", "core") + pct.get_genes("HSF1", "core") + \
        [f"G{i}" for i in range(n_genes - 10)]
    expr = pd.DataFrame(rng.gamma(1.0, 20.0, size=(len(genes), 3)), index=genes,
                        columns=["ct_a", "ct_b", "ct_c"])
    return expr


def test_hpa_reference_scores_pseudobulk_is_mean_log1p():
    expr = _toy_expression()
    out = pct.hpa_reference_scores(expr, scaling="pseudobulk")
    genes = pct.get_genes("NF-κB", "core")
    np.testing.assert_allclose(out["NF-κB"].values, np.log1p(expr.loc[genes]).mean().values)
    assert out.shape == (3, 44)


def test_hpa_reference_scores_single_cell_scaled_normalizes_each_cell_type():
    expr = _toy_expression()
    out = pct.hpa_reference_scores(expr, scaling="single_cell_scaled")
    scaled = expr * (1e4 / expr.sum())
    genes = pct.get_genes("HSF1", "core")
    np.testing.assert_allclose(out["HSF1"].values, np.log1p(scaled.loc[genes]).mean().values)
    # library size matters: doubling one cell type's counts changes nothing after normalization
    expr2 = expr.copy()
    expr2["ct_b"] *= 2
    out2 = pct.hpa_reference_scores(expr2, scaling="single_cell_scaled")
    np.testing.assert_allclose(out["HSF1"].values, out2["HSF1"].values, atol=1e-10)


def test_hpa_reference_scores_warns_on_a_gene_subset():
    expr = _toy_expression(n_genes=400)
    with pytest.warns(UserWarning, match="library_size"):
        pct.hpa_reference_scores(expr, scaling="single_cell_scaled")
    pct.hpa_reference_scores(expr, scaling="single_cell_scaled", library_size=pd.Series(1e6, index=expr.columns))


def test_hpa_reference_scores_rejects_unknown_scaling():
    with pytest.raises(ValueError):
        pct.hpa_reference_scores(_toy_expression(), scaling="bogus")


# ── bundled shift vectors live in the same coordinates ───────────────────────

_REF_DIR = Path(pct.__file__).parent / "reference" / "data"


def _expected(module_deltas):
    return pct.project_shift(pd.Series(module_deltas, dtype=float)).values


def test_attractor_eigenspace_vector_matches_its_module_deltas():
    attr = pct.load_attractor_direction()
    np.testing.assert_allclose(attr["attractor_direction_eigenspace"],
                               _expected(attr["attractor_direction_modules"]), atol=1e-9)


def test_disease_vectors_match_their_module_deltas():
    d = json.loads((_REF_DIR / "disease_vectors.json").read_text(encoding="utf-8"))
    assert d["_meta"]["eigenspace_convention"] == "standardized"
    for name, rec in d.items():
        if name.startswith("_"):
            continue
        np.testing.assert_allclose(rec["mean_vector_12D"],
                                   _expected(rec["module_delta_mean_across_groups"]), atol=1e-9)
        for group in rec["per_cell_type"].values():
            np.testing.assert_allclose(group["eigenspace_vector"], _expected(group["module_deltas"]), atol=1e-9)
            assert group["magnitude"] == pytest.approx(np.linalg.norm(group["eigenspace_vector"]))


def test_aging_vectors_match_their_module_deltas():
    d = json.loads((_REF_DIR / "aging_reference.json").read_text(encoding="utf-8"))
    assert d["eigenspace_convention"] == "standardized"
    for tv in d["tissue_vectors"].values():
        np.testing.assert_allclose(tv["vector_PCs"], _expected(tv["module_deltas"]), atol=1e-9)
    assert d["inflammaging_direction"] == d["tissue_vectors"]["blood"]["vector_PCs"]
    assert d["collapse_direction"] == d["tissue_vectors"]["bone_marrow"]["vector_PCs"]
