"""Eigenspace tests — project + rebuild."""

import numpy as np
import pandas as pd

import perceptome as pct


def test_project_returns_pcs():
    ref = pct.load_hpa_perceptivity()
    sample = ref["R"].iloc[:5]
    out = pct.project(sample)
    assert "coordinates" in out
    assert out["coordinates"].shape[0] == 5
    assert out["coordinates"].shape[1] == 9  # v0.3 = 9 Kaiser PCs


def test_project_confidence_top_pcs_stable():
    ref = pct.load_hpa_perceptivity()
    sample = ref["R"].iloc[:1]
    out = pct.project(sample)
    # PC1-PC6 stable per v0.3 build report
    assert out["confidence"]["PC1"] == "stable"
    assert out["confidence"]["PC2"] == "stable"
    assert out["confidence"]["PC3"] == "stable"


def test_project_handles_missing_modules():
    ref = pct.load_hpa_perceptivity()
    sample = ref["R"].iloc[:3].copy()
    # Drop a column — project should fill with 0
    sample = sample.drop(columns=["NPAS4"])
    out = pct.project(sample)
    assert out["coordinates"].shape == (3, 9)


def test_rebuild_synthetic_matrix():
    rng = np.random.default_rng(0)
    mods = pct.list_modules()
    df = pd.DataFrame(
        rng.normal(size=(60, len(mods))),
        index=[f"ct_{i}" for i in range(60)],
        columns=mods,
    )
    out = pct.rebuild(df, kaiser=True, bootstrap_n=10)
    assert "eigenvalues" in out
    assert "loadings" in out
    assert out["n_modules"] == 44
    assert out["n_cell_types_used"] == 60
