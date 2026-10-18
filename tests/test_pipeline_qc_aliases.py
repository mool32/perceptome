"""Per-cell pipeline (QC -> score -> standardize -> project -> assign), QC rules and gene aliases."""

import warnings

import numpy as np
import pandas as pd
import pytest

import perceptome as pct
from perceptome import convergence as cv


# ── synthetic cells built from the bundled reference profiles ────────────────

def _core_gene_modules():
    cat = pct.load_catalog()["modules"]
    mods_of = {}
    for m, v in cat.items():
        for g in v["core_genes"]:
            mods_of.setdefault(g, []).append(m)
    return sorted(mods_of), mods_of


def make_adata(types, n_per=30, n_filler=2500, seed=0, noise=0.15, rename=None, calibration="pseudobulk"):
    """Cells whose module scores follow the reference profile of ``types`` on the scale of ``calibration``."""
    import anndata
    R = pct.load_calibration(calibration).reference_scores.astype(float)
    genes, mods_of = _core_gene_modules()
    rng = np.random.default_rng(seed)
    rows, labels = [], []
    for t in types:
        target = R.loc[t]
        for _ in range(n_per):
            x = np.array([np.mean([target[m] for m in mods_of[g]]) for g in genes]) + rng.normal(0, noise, len(genes))
            filler = rng.lognormal(0.5, 0.8, n_filler) * (rng.random(n_filler) < 0.7)
            rows.append(np.clip(np.concatenate([x, filler]), 0, None))
            labels.append(t)
    names = [(rename or {}).get(g, g) for g in genes] + [f"FILL{i}" for i in range(n_filler)]
    obs = pd.DataFrame({"truth": labels}, index=[f"c{i}" for i in range(len(rows))])
    return anndata.AnnData(X=np.array(rows, dtype=np.float32), obs=obs, var=pd.DataFrame(index=names))


@pytest.fixture(scope="module")
def cells():
    """Per-cell-scale cells: the pipeline's default calibration is the single-cell-scaled one."""
    return make_adata(["cytotrophoblasts", "megakaryocytes", "hepatocytes", "gastric progenitor cells"],
                      calibration="single_cell_scaled")


def test_assign_cells_recovers_the_cell_type_of_reference_derived_cells(cells):
    res = cv.assign_cells(cells)
    a = res["assignment"]
    truth = cells.obs["truth"].reindex(a.index)
    assert (a["nearest"] == truth).mean() > 0.95
    assert set(res) == {"assignment", "scores", "coordinates", "coverage", "qc", "n_genes"}
    assert (res["n_genes"] >= 1000).all() and list(res["n_genes"].index) == list(a.index)
    assert res["scores"].shape == (len(a), 44) and res["coordinates"].shape == (len(a), 9)
    assert res["qc"]["n_input"] == res["qc"]["n_pass"] == cells.n_obs


def test_core_occupancy_separates_core_from_non_core_cells(cells):
    a = cv.assign_cells(cells)["assignment"]
    truth = cells.obs["truth"].reindex(a.index)
    assert cv.core_fraction(a[truth.isin(["cytotrophoblasts", "megakaryocytes"])]) > 0.95
    assert cv.core_fraction(a[truth == "hepatocytes"]) < 0.05


def test_assign_cells_defaults_to_the_single_cell_calibration_and_the_pseudobulk_one_is_the_wrong_scale(cells):
    """Notes S-B, S-H: per-cell scores against the population-averaged reference are mis-assigned."""
    assert cv.assign_cells(cells)["assignment"].attrs["calibration"] == "single_cell_scaled"
    with pytest.warns(UserWarning, match="single_cell_scaled"):
        wrong = cv.assign_cells(cells, calibration="pseudobulk")["assignment"]
    truth = cells.obs["truth"].reindex(wrong.index)
    assert (wrong["nearest"] == truth).mean() < 0.5


def test_pseudobulk_scale_cells_are_assigned_under_the_pseudobulk_calibration():
    ad = make_adata(["cytotrophoblasts", "hepatocytes"], n_per=20)
    a = cv.assign_cells(ad, calibration="pseudobulk")["assignment"]
    assert (a["nearest"] == ad.obs["truth"].reindex(a.index)).mean() > 0.95


def test_assign_cells_reports_qc_failures():
    ad = make_adata(["hepatocytes"], n_per=20, calibration="single_cell_scaled")
    X = ad.X.copy()
    X[:5, 500:] = 0                       # five cells left with <1,000 detected genes
    ad.X = X
    res = cv.assign_cells(ad)
    assert res["qc"]["n_input"] == 20 and res["qc"]["n_pass"] == 15
    assert len(res["assignment"]) == 15


def test_assign_cells_runs_shape_only(cells):
    a = cv.assign_cells(cells, shape_only=True)["assignment"]
    assert a.attrs["calibration"].endswith("+shape_only")
    assert len(a) == cells.n_obs


# ── QC rules ─────────────────────────────────────────────────────────────────

def test_prepare_expression_applies_the_1000_gene_floor_and_keeps_log_data_untouched():
    ad = make_adata(["hepatocytes"], n_per=10)
    X = ad.X.copy()
    X[0, 300:] = 0
    ad.X = X
    out, rep = cv.prepare_expression(ad)
    assert (rep["n_input"], rep["n_pass"]) == (10, 9) and rep["normalized"] is False
    np.testing.assert_allclose(out.X, ad.X[1:], atol=1e-6)
    assert "n_genes_detected" in out.obs and (out.obs["n_genes_detected"] >= 1000).all()


def test_prepare_expression_normalizes_counts_to_1e4_and_log1p():
    rng = np.random.default_rng(0)
    import anndata
    counts = rng.poisson(2.0, size=(8, 1500)).astype(np.float32)
    counts[:, 0] = 500                      # max > 100 -> treated as counts
    ad = anndata.AnnData(X=counts, var=pd.DataFrame(index=[f"G{i}" for i in range(1500)]))
    out, rep = cv.prepare_expression(ad)
    assert rep["normalized"] is True
    np.testing.assert_allclose(np.expm1(out.X).sum(axis=1), 1e4, rtol=1e-4)


def test_prepare_expression_uppercases_deduplicates_and_keeps_nkx3_1():
    import anndata
    names = ["Nkx3-1", "nkx3-1", "Arntl"] + [f"g{i}" for i in range(1200)]
    ad = anndata.AnnData(X=np.ones((4, len(names)), dtype=np.float32), var=pd.DataFrame(index=names))
    out, _ = cv.prepare_expression(ad)
    assert "NKX3-1" in out.var_names and "NKX3" not in out.var_names     # no '-<digits>' stripping (Note S-Q)
    assert list(out.var_names).count("NKX3-1") == 1 and "ARNTL" in out.var_names


# ── gene aliases (Note S-Q) ──────────────────────────────────────────────────

_RENAMED = {"MB21D1": "CGAS", "TMEM173": "STING1", "ARNTL": "BMAL1"}


def test_alias_table_has_the_three_symbols_the_paper_audited():
    assert {g: pct.aliases_of(g)[0] for g in _RENAMED} == _RENAMED
    assert pct.aliases_of("IRF3") == ()


def test_catalog_still_publishes_the_deprecated_symbols():
    """Table S9 lists MB21D1, TMEM173 and ARNTL; the catalog keeps them verbatim."""
    assert {"MB21D1", "TMEM173"} <= set(pct.get_genes("cGAS-STING", "core"))
    assert "ARNTL" in pct.get_genes("Circadian", "core")


def test_current_symbols_lose_genes_by_exact_matching_and_the_warning_says_so():
    ad = make_adata(["hepatocytes"], n_per=5, rename=_RENAMED)
    with pytest.warns(UserWarning, match="MB21D1 -> CGAS"):
        res = pct.score_modules(ad)
    cov = res["coverage"]
    assert cov.loc["cGAS-STING", "n_genes_found"] == 3 and cov.loc["Circadian", "n_genes_found"] == 4
    assert cov.loc["cGAS-STING", "aliases_used"] == ""


def test_resolve_aliases_recovers_all_five_genes():
    ad = make_adata(["hepatocytes"], n_per=5, rename=_RENAMED)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = pct.score_modules(ad, resolve_aliases=True)
    cov = res["coverage"]
    assert cov.loc["cGAS-STING", "n_genes_found"] == 5 and cov.loc["Circadian", "n_genes_found"] == 5
    assert "MB21D1>CGAS" in cov.loc["cGAS-STING", "aliases_used"]
    assert res["missing_genes"]["cGAS-STING"] == []


def test_either_spelling_scores_the_same_five_genes_with_aliases():
    old = pct.score_modules(make_adata(["oocytes"], n_per=5, noise=0.0), resolve_aliases=True)["scores"]
    new = pct.score_modules(make_adata(["oocytes"], n_per=5, noise=0.0, rename=_RENAMED), resolve_aliases=True)["scores"]
    np.testing.assert_allclose(old.values, new.values, atol=1e-9)


def test_old_symbols_need_no_alias_and_do_not_warn():
    ad = make_adata(["hepatocytes"], n_per=5)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = pct.score_modules(ad)
    assert res["coverage"].loc["cGAS-STING", "n_genes_found"] == 5


def test_resolve_genes_helper():
    found, missing, aliased = pct.resolve_genes(["MB21D1", "IRF3", "NOPE"], {"CGAS", "IRF3"})
    assert found == ["CGAS", "IRF3"] and missing == ["NOPE"] and aliased == {"MB21D1": "CGAS"}
