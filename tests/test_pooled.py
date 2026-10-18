"""Within-patient deepening pooled over datasets with a dataset random intercept (Note S-D)."""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import perceptome as pct
from perceptome import convergence as cv
from perceptome.convergence.pooled import _reml_intercept


def grouped(sizes, mu=0.12, s_u=0.06, s_e=0.05, seed=0):
    rng = np.random.default_rng(seed)
    y, g = [], []
    for i, n in enumerate(sizes):
        y += list(mu + rng.normal(0, s_u) + rng.normal(0, s_e, n))
        g += [f"D{i}"] * n
    return np.array(y), np.array(g)


def test_balanced_design_matches_the_anova_closed_form():
    """Balanced one-way random effects: REML equals the ANOVA estimator while the latter is positive."""
    k, n = 6, 10
    y, g = grouped([n] * k, seed=2)
    groups = y.reshape(k, n)
    msb = n * groups.mean(axis=1).var(ddof=1)
    msw = ((groups - groups.mean(axis=1, keepdims=True)) ** 2).sum() / (k * (n - 1))
    assert msb > msw
    fit = _reml_intercept(y, g)
    assert fit["mu"] == pytest.approx(y.mean(), abs=1e-12)
    assert fit["s2"] == pytest.approx(msw, rel=1e-5)
    assert fit["s_u2"] == pytest.approx((msb - msw) / n, rel=1e-4)
    assert fit["se"] == pytest.approx(np.sqrt(msb / (k * n)), rel=1e-5)


def test_without_between_dataset_variance_it_reduces_to_the_ordinary_standard_error():
    rng = np.random.default_rng(1)
    parts = []
    for n in (12, 20, 9, 15):
        e = rng.normal(0, 0.05, n)
        parts.append(0.1 + e - e.mean())             # identical group means -> no between-group variance
    y, g = np.concatenate(parts), np.repeat(["a", "b", "c", "d"], (12, 20, 9, 15))
    fit = _reml_intercept(y, g)
    assert fit["s_u2"] == 0.0
    assert fit["mu"] == pytest.approx(0.1)
    assert fit["se"] == pytest.approx(y.std(ddof=1) / np.sqrt(len(y)), rel=1e-9)


def test_one_dataset_is_the_plain_mean():
    d = pd.Series(np.random.default_rng(0).normal(0.1, 0.05, 25))
    r = cv.pooled_within_patient({"only": d})
    assert r["n_datasets"] == 1 and r["n_patients"] == 25
    assert r["estimate"] == pytest.approx(d.mean())
    assert r["se"] == pytest.approx(d.std(ddof=1) / 5)


def test_the_mixed_model_credits_between_dataset_clustering():
    """Heterogeneous datasets: the pooled standard error is larger than the one that ignores the clustering."""
    y, g = grouped([30, 12, 25, 8], mu=0.1, s_u=0.15, s_e=0.04, seed=5)
    fit = _reml_intercept(y, g)
    assert fit["s_u2"] > 0
    assert fit["se"] > y.std(ddof=1) / np.sqrt(len(y)) * 1.5


@pytest.mark.filterwarnings("ignore:The MLE may be on the boundary")
def test_it_agrees_with_statsmodels_where_the_solution_is_interior():
    smf = pytest.importorskip("statsmodels.formula.api")
    y, g = grouped([14, 30, 22, 9, 17, 26, 11, 19], mu=0.12, s_u=0.08, s_e=0.06, seed=11)
    fit = _reml_intercept(y, g)
    df = pd.DataFrame({"y": y, "g": g})
    m = smf.mixedlm("y ~ 1", df, groups=df["g"]).fit(reml=True)
    assert m.converged and fit["s_u2"] > 0
    assert fit["mu"] == pytest.approx(m.fe_params["Intercept"], abs=1e-7)
    assert fit["se"] == pytest.approx(m.bse_fe["Intercept"], rel=2e-3)
    assert fit["s2"] == pytest.approx(m.scale, rel=1e-4)
    assert fit["s_u2"] == pytest.approx(float(m.cov_re.iloc[0, 0]), rel=2e-3)


def test_input_forms_agree_and_the_one_sided_p_is_the_papers_convention():
    y, g = grouped([20, 25, 8, 12], seed=3)
    long = pd.DataFrame({"delta": y, "dataset": g})
    as_dict = {name: sub[["delta"]] for name, sub in long.groupby("dataset")}
    a = cv.pooled_within_patient(long, dataset="dataset")
    b = cv.pooled_within_patient(as_dict)
    assert a["estimate"] == pytest.approx(b["estimate"]) and a["se"] == pytest.approx(b["se"])
    assert a["estimate_points"] == pytest.approx(100 * a["estimate"])
    assert a["z"] == pytest.approx(a["estimate"] / a["se"])
    assert a["p_one_sided"] == pytest.approx(stats.norm.sf(a["z"]))
    assert a["p_two_sided"] == pytest.approx(2 * stats.norm.sf(abs(a["z"])))
    assert pct.pooled_within_patient is cv.pooled_within_patient
    # Note S-D: +13.0 points, SE 1.87 -> z = 6.95, one-sided P = 1.9e-12
    assert 13.0 / 1.87 == pytest.approx(6.95, abs=0.005)
    assert stats.norm.sf(6.95) == pytest.approx(1.9e-12, rel=0.1)


def test_it_takes_the_per_patient_frames_of_within_patient_occupancy_and_summarizes_datasets():
    def frame(hits_m, hits_n, n=40):
        labels, patient, arm = [], [], []
        for p, (hm, hn) in enumerate(zip(hits_m, hits_n)):
            for a, hits in (("tumour", hm), ("normal", hn)):
                labels += ["cytotrophoblasts"] * hits + ["hepatocytes"] * (n - hits)
                patient += [f"P{p}"] * n
                arm += [a] * n
        idx = pd.RangeIndex(len(labels))
        return pd.Series(labels, idx), pd.Series(patient, idx), pd.Series(arm, idx)

    runs = {"colon": ([16, 14, 18, 12, 20, 15], [4, 3, 5, 4, 2, 4]),
            "liver": ([10, 9, 3], [4, 4, 6])}
    frames = {}
    for name, (m, n) in runs.items():
        near, pat, arm = frame(m, n)
        frames[name] = cv.within_patient_occupancy(near, pat, arm, "tumour", "normal", region="anchor")["per_patient"]
    r = cv.pooled_within_patient(frames)
    assert r["n_patients"] == 9 and r["n_datasets"] == 2
    assert list(r["per_dataset"]["n_patients"]) == [6, 3] and list(r["per_dataset"]["n_deeper"]) == [6, 2]
    assert r["estimate_points"] > 0
    matched = cv.pooled_within_patient({k: v.assign(delta_matched=v["delta"]) for k, v in frames.items()},
                                       column="delta_matched")
    assert matched["estimate"] == pytest.approx(r["estimate"])


def test_input_validation():
    with pytest.raises(ValueError, match="dataset"):
        cv.pooled_within_patient(pd.DataFrame({"delta": [0.1, 0.2, 0.3]}))
    with pytest.raises(ValueError, match="at least two"):
        cv.pooled_within_patient({"a": pd.Series([0.1])})
    with pytest.raises(ValueError, match="no column"):
        cv.pooled_within_patient({"a": pd.DataFrame({"other": [0.1, 0.2]})})
