"""Perceptivity tests — HPA reference, capacity floor, predict_engagement."""

import numpy as np
import pandas as pd
import pytest

import perceptome as pct


def test_hpa_reference_dimensions():
    ref = pct.load_hpa_perceptivity()
    assert ref["R"].shape == (154, 44)
    assert ref["A"].shape == (154, 44)
    assert ref["C"].shape == (154, 44)
    assert ref["headroom"].shape == (154, 44)
    assert "NPAS4" in ref["R"].columns
    assert "brain excitatory neurons" in ref["R"].index


def test_npas4_neuron_enriched():
    """NPAS4 A_baseline should be highest in neuronal cell type class."""
    ref = pct.load_hpa_perceptivity()
    A_npas4 = ref["A"]["NPAS4"]
    klass = ref["cell_type_class"]
    by_class = pd.DataFrame({"A": A_npas4, "klass": klass}).groupby("klass")["A"].mean()
    # Neurons should be top or near-top
    top = by_class.sort_values(ascending=False).head(2).index
    assert any("neuron" in t.lower() for t in top)


def test_capacity_floor_classifications():
    assert pct.capacity_floor(5.5) == "saturated_blocked_up"
    assert pct.capacity_floor(2.0) == "capacious"
    assert pct.capacity_floor(3.5) == "intermediate"
    assert pct.capacity_floor(float("nan")) == "no_data"
    # boundary cases: lo=2.5, hi=4.5
    assert pct.capacity_floor(4.5) == "intermediate"   # > 4.5 saturated, == intermediate
    assert pct.capacity_floor(4.51) == "saturated_blocked_up"
    assert pct.capacity_floor(2.49) == "capacious"


def test_predict_engagement_paper45_goblet():
    """Paper 4.5 a-priori: enteric stem cells UPR-ATF6 should be saturated."""
    pred = pct.predict_engagement("enteric stem cells", ["UPR-ATF6"])
    assert pred.loc["UPR-ATF6", "capacity_floor"] == "saturated_blocked_up"
    assert pred.loc["UPR-ATF6", "A_baseline"] > 4.5


def test_predict_engagement_cardiomyocyte_hsf1():
    """Paper 4.4 a-priori: cardiomyocyte HSF1 should be capacious."""
    pred = pct.predict_engagement("cardiomyocytes", ["HSF1"])
    assert pred.loc["HSF1", "capacity_floor"] == "capacious"
    assert pred.loc["HSF1", "A_baseline"] < 2.5


def test_predict_engagement_unknown_cell_type_raises():
    with pytest.raises(KeyError):
        pct.predict_engagement("xyz_nonexistent_cell_type", ["HSF1"])


def test_compute_perceptivity_with_synthetic():
    ref = pct.load_hpa_perceptivity()
    R = ref["R"].iloc[:10]
    A = ref["A"].iloc[:10]
    out = pct.compute_perceptivity(R, A)
    assert "per_cell_type" in out
    assert "per_module" in out
    assert out["per_module"]["C"].shape == (10, 44)
    assert "spec_quadrant" in out["per_cell_type"].columns


def test_quadrant_classifier():
    bs_med, gs_med = 5.0, 5.0
    assert "Q1" in pct.classify_quadrant(7, 7, bs_med, gs_med)
    assert "Q2" in pct.classify_quadrant(3, 7, bs_med, gs_med)
    assert "Q3" in pct.classify_quadrant(7, 3, bs_med, gs_med)
    assert "Q4" in pct.classify_quadrant(3, 3, bs_med, gs_med)
    assert pct.classify_quadrant(float("nan"), 5, bs_med, gs_med) == "NA"
