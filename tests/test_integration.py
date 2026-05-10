"""End-to-end integration tests + reproducibility regressions."""

import numpy as np
import pandas as pd

import perceptome as pct


def test_pipeline_score_project_compare(tiny_adata):
    """Score → project → compare_to_references runs end-to-end without error."""
    scores = pct.score_modules(tiny_adata, method="mean_raw")["scores"]
    coords = pct.project(scores)["coordinates"]
    refs = pct.compare_to_references(coords, include_diseases=False, include_aging=False)
    assert "cos_attractor" in refs
    assert "eigenspace_region" in refs
    assert len(refs["cos_attractor"]) == tiny_adata.n_obs


def test_pipeline_perceptivity_predict():
    """Load HPA → predict_engagement returns per-module capacity prediction."""
    pred = pct.predict_engagement(
        "hepatocytes",
        ["UPR-ATF6", "HSF1", "Cell Cycle"],
    )
    assert set(pred.index) == {"UPR-ATF6", "HSF1", "Cell Cycle"}
    for col in ("R_baseline", "A_baseline", "C", "headroom",
                "capacity_floor", "predicted_direction"):
        assert col in pred.columns


def test_reproducibility_score_modules(tiny_adata):
    """Same input → same output, deterministic across two runs."""
    s1 = pct.score_modules(tiny_adata, method="mean_raw")["scores"]
    s2 = pct.score_modules(tiny_adata, method="mean_raw")["scores"]
    pd.testing.assert_frame_equal(s1, s2)


def test_reproducibility_eigenspace_projection():
    """Project the same scores twice → identical coordinates."""
    ref = pct.load_hpa_perceptivity()
    sample = ref["R"].iloc[:5]
    c1 = pct.project(sample)["coordinates"]
    c2 = pct.project(sample)["coordinates"]
    pd.testing.assert_frame_equal(c1, c2)


def test_reproducibility_perceptivity():
    """compute_perceptivity is deterministic for the same R/A inputs."""
    ref = pct.load_hpa_perceptivity()
    R, A = ref["R"].iloc[:20], ref["A"].iloc[:20]
    p1 = pct.compute_perceptivity(R, A)
    p2 = pct.compute_perceptivity(R, A)
    pd.testing.assert_frame_equal(p1["per_cell_type"], p2["per_cell_type"])
    pd.testing.assert_frame_equal(p1["per_module"]["C"], p2["per_module"]["C"])


def test_hpa_reference_signs_consistent_with_v1_4():
    """Spot-check three v1.4-locked numbers in the HPA reference.

    These are from the closed metric design v1.4 (perceptivity_metric/results/).
    If they drift, something changed about the underlying HPA matrix or scoring.
    """
    ref = pct.load_hpa_perceptivity()
    A = ref["A"]
    # Plasmablast precursor: 'b-cells' (HPA naming, not 'naive b-cells')
    assert 4.0 < A.loc["b-cells", "UPR-ATF6"] < 5.0  # observed 4.47
    # Hepatocytes saturated UPR-ATF6
    assert A.loc["hepatocytes", "UPR-ATF6"] > 6.0  # observed 6.32
    # Cardiomyocytes capacious HSF1
    assert A.loc["cardiomyocytes", "HSF1"] < 2.5  # observed 2.03


def test_npas4_zero_in_non_neuron_synthetic():
    """Synthetic AnnData has neuron-specific genes — NPAS4 score should be measurable
    even when filled uniformly; this is a regression check that the catalog gene list
    actually maps."""
    import anndata
    rng = np.random.default_rng(0)
    genes = pct.get_genes("NPAS4", "core") + pct.get_genes("NPAS4", "activity")
    X = rng.normal(loc=2.0, scale=0.5, size=(20, len(genes)))
    a = anndata.AnnData(
        X=np.clip(X, 0, None),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(20)]),
        var=pd.DataFrame(index=genes),
    )
    res = pct.score_modules(a, gene_set="core", method="mean_raw")
    assert res["coverage"].loc["NPAS4", "n_genes_found"] == 5
    assert (res["scores"]["NPAS4"] > 0).any()


def test_attractor_alignment_consistent_in_both_modes():
    """Eigenspace and modules modes both use the same underlying direction; aligned
    inputs should pass the P3 threshold in both modes."""
    attr = pct.load_attractor_direction()
    out_e = pct.attractor_alignment(attr["attractor_direction_eigenspace"], mode="eigenspace")
    out_m = pct.attractor_alignment(attr["attractor_direction_modules"], mode="modules")
    assert out_e["passes_p3_threshold"]
    assert out_m["passes_p3_threshold"]


def test_paper4_5_paper4_4_anchors_predict_correctly():
    """Bundled regression: 4 a-priori predictions that paper4.5 / paper4.4 reported."""
    expected = [
        # (cell_type, module, expected_floor)
        ("enteric stem cells", "UPR-ATF6", "saturated_blocked_up"),  # paper4.5 PASS
        ("hepatocytes", "UPR-ATF6", "saturated_blocked_up"),         # paper4.7 cont.
        ("cardiomyocytes", "HSF1", "capacious"),                      # paper4.4 PASS
        ("brain excitatory neurons", "HSF1", "capacious"),            # paper4 N6
    ]
    for ct, mod, expected_floor in expected:
        pred = pct.predict_engagement(ct, [mod])
        assert pred.loc[mod, "capacity_floor"] == expected_floor, (
            f"{ct}/{mod}: expected {expected_floor}, got {pred.loc[mod, 'capacity_floor']} "
            f"(A_baseline={pred.loc[mod, 'A_baseline']:.2f})"
        )
