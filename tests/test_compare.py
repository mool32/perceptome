"""compare/ tests — conditions, regime, references, divergence."""

import numpy as np
import pandas as pd
import pytest

import perceptome as pct


# ── infrastructure_regime ─────────────────────────────────────────────────────

def _delta_series(perc_delta, infra_delta):
    """Build a per-module delta Series with given perception/infra means."""
    perc = ["NF-κB", "ERK/MAPK", "JAK-STAT", "cAMP/CREB", "NFAT", "Wnt"]
    infra = ["HSF1", "UPR-IRE1", "UPR-PERK", "UPR-ATF6", "NRF2", "mTOR"]
    return pd.Series({**{m: perc_delta for m in perc}, **{m: infra_delta for m in infra}})


def test_infrastructure_regime_supply_chain():
    out = pct.infrastructure_regime(_delta_series(+0.5, +0.5))
    assert out["regime"] == "supply_chain"


def test_infrastructure_regime_firefighting():
    out = pct.infrastructure_regime(_delta_series(-0.5, +0.5))
    assert out["regime"] == "firefighting"


def test_infrastructure_regime_collapse():
    out = pct.infrastructure_regime(_delta_series(-0.5, -0.5))
    assert out["regime"] == "collapse"


def test_infrastructure_regime_unsupported():
    out = pct.infrastructure_regime(_delta_series(+0.5, -0.5))
    assert out["regime"] == "unsupported"


def test_infrastructure_regime_dataframe_input():
    df = pd.DataFrame({"delta": _delta_series(+0.5, +0.5)})
    out = pct.infrastructure_regime(df)
    assert out["regime"] == "supply_chain"
    assert out["perception_delta"] > 0
    assert out["infrastructure_delta"] > 0


def test_infrastructure_regime_includes_npas4():
    """NPAS4 is in the v0.2 perception list — its delta should affect perception_delta."""
    s = pd.Series({"NPAS4": +1.0, "HSF1": 0.0})
    out = pct.infrastructure_regime(s)
    assert out["perception_delta"] > 0


# ── compare_conditions ────────────────────────────────────────────────────────

def _make_two_conditions():
    """Two synthetic score DataFrames: cond2 has +1.0 shift on first 10 modules."""
    rng = np.random.default_rng(0)
    modules = pct.list_modules()
    n = 100
    s1 = pd.DataFrame(rng.normal(size=(n, len(modules))), columns=modules,
                      index=[f"c1_{i}" for i in range(n)])
    s2 = s1.copy()
    s2.index = [f"c2_{i}" for i in range(n)]
    s2.iloc[:, :10] += 1.0  # shift first 10 modules up
    return s1, s2


def test_compare_conditions_returns_expected_keys():
    s1, s2 = _make_two_conditions()
    out = pct.compare_conditions(s2, s1, project_to_eigenspace=True)
    for k in ("delta_modules", "delta_eigenspace", "magnitude",
              "top_modules", "disease_type", "infrastructure_regime"):
        assert k in out


def test_compare_conditions_detects_shift():
    """First 10 modules shifted +1.0; mean delta on those should be positive."""
    s1, s2 = _make_two_conditions()
    out = pct.compare_conditions(s2, s1, project_to_eigenspace=False)
    df = out["delta_modules"]
    shifted_mods = list(s1.columns[:10])
    assert df.loc[shifted_mods, "delta"].mean() > 0.5


def test_compare_conditions_top_modules_picks_largest():
    s1, s2 = _make_two_conditions()
    out = pct.compare_conditions(s2, s1, project_to_eigenspace=False)
    top = out["top_modules"]
    assert len(top) == 10
    # the shifted modules should dominate top-10 by |delta|
    shifted = set(s1.columns[:10])
    assert len(set(top.index) & shifted) >= 5


def test_compare_conditions_per_cell_type():
    s1, s2 = _make_two_conditions()
    obs1 = pd.Series(["A"] * 50 + ["B"] * 50, index=s1.index)
    obs2 = pd.Series(["A"] * 50 + ["B"] * 50, index=s2.index)
    out = pct.compare_conditions(
        s2, s1, cell_type_column="cell_type", obs1=obs2, obs2=obs1,
        project_to_eigenspace=True,
    )
    assert "per_cell_type" in out
    assert "A" in out["per_cell_type"] and "B" in out["per_cell_type"]
    assert "_cross_cell_type" in out["per_cell_type"]


# ── divergence_score ──────────────────────────────────────────────────────────

def test_divergence_score_classifies_convergent():
    """Manually build per_cell_type dict with same eigenspace_vector for immune+structural."""
    per_ct = {
        "macrophages": {"eigenspace_vector": [1.0, 0.5, 0.0]},
        "fibroblasts": {"eigenspace_vector": [1.0, 0.5, 0.0]},
    }
    out = pct.divergence_score(per_ct)
    assert out["pattern"] == "convergent"
    assert out["divergence_score"] > 0.99


def test_divergence_score_classifies_divergent():
    per_ct = {
        "macrophages": {"eigenspace_vector": [1.0, 0.5, 0.0]},
        "fibroblasts": {"eigenspace_vector": [-1.0, -0.5, 0.0]},
    }
    out = pct.divergence_score(per_ct)
    assert out["pattern"] == "divergent"
    assert out["divergence_score"] < -0.99


def test_divergence_score_insufficient_data():
    per_ct = {"macrophages": {"eigenspace_vector": [1.0, 0.0, 0.0]}}
    out = pct.divergence_score(per_ct)
    assert out["pattern"] == "insufficient_data"


# ── compare_to_references ─────────────────────────────────────────────────────

def test_compare_to_references_attractor_cosine():
    """Coordinates aligned with attractor direction should give high cos_attractor."""
    attr = pct.load_attractor_direction()
    direction = attr["attractor_direction_eigenspace"]
    coords = pd.DataFrame([direction, -direction, np.zeros_like(direction)],
                          index=["aligned", "anti", "zero"],
                          columns=[f"PC{i+1}" for i in range(len(direction))])
    out = pct.compare_to_references(coords, include_diseases=False, include_aging=False)
    assert "cos_attractor" in out
    assert out["cos_attractor"]["aligned"] > 0.99
    assert out["cos_attractor"]["anti"] < -0.99


def test_compare_to_references_disease_vectors_recomputed():
    """Disease vectors recomputed on v0.3 eigenspace — should return DataFrame of cosines."""
    attr = pct.load_attractor_direction()
    direction = attr["attractor_direction_eigenspace"]
    coords = pd.DataFrame([direction], index=["sample"],
                          columns=[f"PC{i+1}" for i in range(len(direction))])
    out = pct.compare_to_references(coords, include_diseases=True, include_aging=False,
                                    include_attractor=False)
    assert "disease_similarities" in out
    df = out["disease_similarities"]
    # Should have RA/AD/IPF/DKD as columns (4 diseases recomputed)
    expected = {"RA", "AD", "IPF", "DKD"}
    assert expected.issubset(set(df.columns)), f"Missing diseases: {expected - set(df.columns)}"


def test_compare_to_references_aging_recomputed():
    """Aging vectors recomputed — should return inflammaging + collapse cosines."""
    attr = pct.load_attractor_direction()
    direction = attr["attractor_direction_eigenspace"]
    coords = pd.DataFrame([direction], index=["sample"],
                          columns=[f"PC{i+1}" for i in range(len(direction))])
    out = pct.compare_to_references(coords, include_diseases=False, include_aging=True,
                                    include_attractor=False)
    assert "cos_inflammaging" in out
    assert "cos_collapse" in out
    assert "aging_type" in out


def test_eigenspace_region_classifier():
    """compare_to_references always emits an eigenspace_region label per sample."""
    coords = pd.DataFrame({
        "PC1": [3.0, -2.0, 0.0, 0.0],
        "PC2": [0.0, 0.0, 2.0, -2.0],
    }, index=["high_perc", "low_perc", "stress", "growth"])
    out = pct.compare_to_references(coords, include_diseases=False, include_aging=False,
                                    include_attractor=False)
    region = out["eigenspace_region"]
    assert region["high_perc"] == "high_perception"
    assert region["low_perc"] == "low_perception"
    assert region["stress"] == "stress_response"
    assert region["growth"] == "growth_mode"
