"""Attractor reference tests — Paper 4.2 P3 capacity-direction."""

import numpy as np
import pandas as pd

import perceptome as pct


def test_load_attractor_direction():
    attr = pct.load_attractor_direction()
    assert "attractor_cluster_cells" in attr
    assert len(attr["attractor_cluster_cells"]) == 8
    assert "gastric chief cells" in attr["attractor_cluster_cells"]
    assert isinstance(attr["attractor_direction_modules"], pd.Series)
    assert isinstance(attr["attractor_direction_eigenspace"], np.ndarray)


def test_attractor_direction_signs():
    """Top up modules per Paper 4.2 closing memo: Cell Cycle, Hippo, UPR-ATF6, HIF.
    Top down: GR, MR, AR, cAMP/CREB, LXR.
    """
    attr = pct.load_attractor_direction()
    deltas = attr["attractor_direction_modules"]
    for up_mod in ("Cell Cycle", "Hippo", "UPR-ATF6", "HIF"):
        assert deltas[up_mod] > 0, f"{up_mod} expected +; got {deltas[up_mod]}"
    for down_mod in ("GR", "MR", "AR", "cAMP/CREB", "LXR"):
        assert deltas[down_mod] < 0, f"{down_mod} expected -; got {deltas[down_mod]}"


def test_attractor_alignment_module_mode():
    """A shift vector identical to attractor direction should give cosine ≈ 1."""
    attr = pct.load_attractor_direction()
    direction = attr["attractor_direction_modules"]
    out = pct.attractor_alignment(direction, mode="modules")
    assert out["cosine"] > 0.99
    assert out["passes_p3_threshold"] is True


def test_attractor_alignment_eigenspace_mode():
    """A shift vector identical to attractor eigenspace should give cosine ≈ 1."""
    attr = pct.load_attractor_direction()
    out = pct.attractor_alignment(attr["attractor_direction_eigenspace"], mode="eigenspace")
    assert out["cosine"] > 0.99


def test_attractor_alignment_orthogonal_low_cosine():
    """A shift vector opposite to attractor should give cosine ≈ -1."""
    attr = pct.load_attractor_direction()
    direction = attr["attractor_direction_modules"]
    out = pct.attractor_alignment(-direction, mode="modules")
    assert out["cosine"] < -0.99
    assert out["passes_p3_threshold"] is False
