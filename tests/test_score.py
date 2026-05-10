"""Scoring tests — mean_raw / mean_zscore / scanpy_corrected."""

import numpy as np
import pandas as pd

import perceptome as pct


def test_score_modules_returns_44_columns(tiny_adata):
    res = pct.score_modules(tiny_adata, method="mean_raw")
    assert res["scores"].shape == (50, 44)
    assert "NPAS4" in res["scores"].columns


def test_score_modules_has_expected_keys(tiny_adata):
    res = pct.score_modules(tiny_adata)
    for k in ("scores", "coverage", "missing_genes", "cognitive_load"):
        assert k in res
    assert res["coverage"].shape[0] == 44


def test_mean_raw_preserves_absolute(tiny_adata):
    res = pct.score_modules(tiny_adata, method="mean_raw")
    # mean expression of NPAS4 core genes should be positive (we filled with positives)
    assert res["scores"]["NPAS4"].mean() > 0


def test_mean_zscore_normalizes(tiny_adata):
    res = pct.score_modules(tiny_adata, method="mean_zscore")
    # for any module with non-zero std, mean across cells ≈ 0
    for mod in res["scores"].columns:
        col = res["scores"][mod]
        if col.std() > 1e-6:
            assert abs(col.mean()) < 1e-9


def test_score_readiness_vs_activity_differ(tiny_adata):
    R = pct.score_readiness(tiny_adata)
    A = pct.score_activity(tiny_adata)
    assert R.shape == A.shape == (50, 44)
    # for at least one module they must differ (different gene sets)
    assert (R - A).abs().sum().sum() > 0


def test_missing_module_genes_filled_with_zero(tiny_adata):
    res = pct.score_modules(tiny_adata, gene_set="core", min_genes=2)
    # Modules without enough gene coverage in the tiny AnnData are filled with 0
    cov = res["coverage"]
    low_cov = cov[cov["n_genes_found"] < 2]
    for mod in low_cov.index:
        assert (res["scores"][mod] == 0).all()
