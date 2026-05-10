"""Drug subpackage tests — anchors + activity_layer_screen."""

import numpy as np
import pandas as pd
import pytest

import perceptome as pct


def test_drug_anchors_returns_9():
    df = pct.drug_anchors()
    assert len(df) == 9
    rescues = pct.drug_anchors(role="validated_rescue")
    controls = pct.drug_anchors(role="positive_control")
    assert len(rescues) == 6
    assert len(controls) == 3


def test_drug_anchors_mek_pair():
    df = pct.drug_anchors()
    mek = df[(df["class"] == "MEKi") & (df["module"] == "ERK/MAPK")]
    assert len(mek) == 1
    assert mek.iloc[0]["expected_sign"] == -1
    assert mek.iloc[0]["block5_z"] == -1.026
    assert mek.iloc[0]["block2_holdout_pass"] is True


def test_drug_anchors_role_filter_invalid_returns_empty():
    df = pct.drug_anchors(role="invalid")
    assert len(df) == 0


def test_activity_layer_screen_runs(tiny_adata):
    """Smoke test — small synthetic dataset, expect graceful behavior."""
    # Need a larger synthetic dataset for the bg pool ≥ 30 requirement
    import anndata
    rng = np.random.default_rng(0)
    n = 100
    genes = ["FOS", "EGR1", "DUSP1", "DUSP6", "MYC", "HSPA1A", "HSPA1B",
             "ACTB", "GAPDH"] + [f"G{i}" for i in range(20)]
    X = rng.normal(loc=0, scale=1, size=(n, len(genes)))
    drugs = ["trametinib"] * 5 + ["selumetinib"] * 5 + [f"drug_{i}" for i in range(40)]
    drugs = drugs + drugs  # 100 rows
    obs = pd.DataFrame({"pert": drugs[:n]}, index=[f"sig_{i}" for i in range(n)])
    a = anndata.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))

    out = pct.activity_layer_screen(
        a, pert_col="pert",
        test_perturbations={"MEKi_test": ["trametinib", "selumetinib"]},
        panels=[("MEKi", "ERK/MAPK", -1)],
        n_perm=200, seed=42,
    )
    assert "observed_z" in out.columns
    assert "p_one_sided" in out.columns
    assert "verdict" in out.columns
    assert len(out) == 1


def test_activity_layer_screen_too_small_bg_raises():
    import anndata
    rng = np.random.default_rng(0)
    a = anndata.AnnData(
        X=rng.normal(size=(10, 5)),
        obs=pd.DataFrame({"pert": ["drugA"] * 5 + ["drugB"] * 5},
                         index=[f"s{i}" for i in range(10)]),
        var=pd.DataFrame(index=["FOS", "EGR1", "DUSP1", "MYC", "ACTB"]),
    )
    with pytest.raises(ValueError, match="Background pool too small"):
        pct.activity_layer_screen(a, pert_col="pert",
                                   test_perturbations="drugA",
                                   panels=[("test", "ERK/MAPK", -1)])


def test_compare_to_references_no_drug_kwarg():
    """include_drugs was removed in v0.2 — passing it should raise TypeError."""
    import inspect
    sig = inspect.signature(pct.compare_to_references)
    assert "include_drugs" not in sig.parameters
