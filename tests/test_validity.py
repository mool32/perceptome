"""Validity scorecard tests."""

import perceptome as pct


def test_random_panel_reproducible(tiny_adata):
    p1 = pct.random_200_panel(tiny_adata, n=20, seed=42)
    p2 = pct.random_200_panel(tiny_adata, n=20, seed=42)
    assert p1 == p2


def test_random_panel_changes_with_seed(tiny_adata):
    p1 = pct.random_200_panel(tiny_adata, n=20, seed=42)
    p2 = pct.random_200_panel(tiny_adata, n=20, seed=999)
    assert p1 != p2


def test_housekeeping_panel_loads():
    hk = pct.housekeeping_panel()
    assert "ACTB" in hk
    assert "GAPDH" in hk


def test_cell_cycle_panel_loads():
    cc = pct.cell_cycle_panel()
    assert "MKI67" in cc
    assert "CCNB1" in cc


def test_validate_perturbation_runs(tiny_adata):
    sc = pct.validate_perturbation(
        tiny_adata, condition_col="condition",
        perturbation_value="drug", control_value="ctrl",
        require_cell_cycle=False,
    )
    assert sc.overall_verdict in ("PASS", "MIXED", "ARTIFACT", "INCONCLUSIVE")
    names = [c.name for c in sc.checks]
    assert "random_200" in names
    assert "housekeeping" in names
    assert "cell_cycle" in names


def test_scorecard_pretty_prints(tiny_adata):
    sc = pct.validate_perturbation(
        tiny_adata, condition_col="condition",
        perturbation_value="drug", control_value="ctrl",
    )
    out = pct.scorecard(sc)
    assert "Overall:" in out
    assert "random_200" in out
