"""Malignant versus matched-normal occupants of the state, module by module (Fig. 6, Fig. S18)."""

import numpy as np
import pandas as pd
import pytest

import perceptome as pct
from perceptome import convergence as cv
from perceptome.convergence.contrast import _bh

MODULES = pct.list_modules()


def occupants(n_patients=8, n_occ=40, n_other=15, seed=0, effect=None, depth_module=None, signature=None,
              short=()):
    """Synthetic cells: module scores around 0, ``n_occ`` core occupants and ``n_other`` others per arm.

    ``effect``: {module: shift} added to malignant occupants. ``depth_module``: a module that rises with the
    number of detected genes in every cell (a depth artifact). ``signature``: a module vector added (x0.5)
    to malignant occupants. ``short``: patients that get only 5 normal occupants.
    """
    rng = np.random.default_rng(seed)
    rows, meta, nearest, genes = [], [], [], []
    for p in range(n_patients):
        for arm in ("normal", "tumour"):
            n_here = 5 if (arm == "normal" and p in short) else n_occ
            for occ in [True] * n_here + [False] * n_other:
                g = int(rng.integers(2000, 6000) if arm == "tumour" else rng.integers(2000, 3500))
                s = pd.Series(rng.normal(0, 0.1, len(MODULES)), index=MODULES)
                if occ and arm == "tumour":
                    for m, v in (effect or {}).items():
                        s[m] += v
                    if signature is not None:
                        s += 0.5 * signature.reindex(MODULES).values
                if depth_module:
                    s[depth_module] += 0.0004 * (g - 3000)
                rows.append(s.values)
                genes.append(g)
                meta.append((f"P{p:02d}", arm))
                nearest.append("cytotrophoblasts" if occ else "hepatocytes")
    idx = [f"c{i}" for i in range(len(rows))]
    info = pd.DataFrame(meta, index=idx, columns=["patient", "arm"])
    return (pd.DataFrame(rows, index=idx, columns=MODULES), pd.Series(nearest, index=idx),
            info["patient"], info["arm"], pd.Series(genes, index=idx))


def run(data, **kw):
    scores, near, patient, arm, genes = data
    kw.setdefault("n_genes", genes)
    return cv.occupant_contrast(scores, near, patient, arm, "tumour", "normal", **kw)


def test_a_module_that_separates_the_occupants_leads_the_table():
    r = run(occupants(effect={"Cell Cycle": 0.5, "mTOR": 0.3}), balance_depth=False)
    t = r["table"]
    assert set(t.index[:2]) == {"Cell Cycle", "mTOR"}
    assert t.loc["Cell Cycle", "n_higher"] == 8 and t.loc["Cell Cycle", "delta"] == pytest.approx(0.5, abs=0.05)
    assert t.loc["Cell Cycle", "wilcoxon_p"] == pytest.approx(2 / 2 ** 8)           # 8 of 8 patients: 2/256
    rest = t.drop(["Cell Cycle", "mTOR"])
    assert rest["wilcoxon_p"].median() > 0.2
    assert (t["q_bh"] >= t["wilcoxon_p"]).all() and t["q_bh"].max() <= 1.0
    assert r["n_patients"] == 8 and list(r["n_cells"].columns) == ["n_malignant", "n_normal"]
    assert r["per_patient"].shape == (8, len(MODULES))


def test_the_depth_balance_removes_a_difference_that_only_follows_library_size():
    data = occupants(depth_module="NRF2", effect={"Cell Cycle": 0.5}, seed=3)
    raw = run(data, balance_depth=False)
    bal = run(data, balance_depth=True)
    assert raw["table"].loc["NRF2", "delta"] > 0.3 and raw["table"].loc["NRF2", "n_higher"] == 8
    assert abs(bal["table"].loc["NRF2", "delta"]) < 0.1
    assert bal["table"].loc["NRF2", "delta_raw"] == pytest.approx(raw["table"].loc["NRF2", "delta_raw"])
    assert bal["table"].loc["Cell Cycle", "delta"] == pytest.approx(0.5, abs=0.06)       # a real difference survives
    d = bal["depth"]
    assert d["malignant_raw_median"] > d["normal_raw_median"] + 800
    assert abs(d["malignant_balanced_median"] - d["normal_balanced_median"]) < 0.1 * d["normal_balanced_median"]
    assert raw["depth"].keys() == {"malignant_raw_median", "normal_raw_median"}


def test_the_difference_correlates_with_the_core_signature_only_if_it_lies_along_it():
    sig = cv.core_signature("pseudobulk")
    along = run(occupants(signature=sig, seed=1), balance_depth=False)
    against = run(occupants(signature=-sig, seed=1), balance_depth=False)
    off = run(occupants(effect={"Cell Cycle": 0.5, "HIF": 0.4}, seed=1), balance_depth=False)
    assert along["rho_core_signature"] > 0.9 and against["rho_core_signature"] < -0.9
    assert abs(off["rho_core_signature"]) < 0.6                  # a handful of modules off the defining axis


def test_patients_below_the_cell_floor_are_left_out():
    r = run(occupants(short=(1, 4)), balance_depth=False)
    assert r["n_patients"] == 6 and "P01" not in r["patients"] and "P04" not in r["patients"]
    low = run(occupants(short=(1, 4)), balance_depth=False, min_cells=5)
    assert low["n_patients"] == 8


def test_median_statistic_and_validation():
    data = occupants(effect={"Cell Cycle": 0.5}, seed=2)
    r = run(data, balance_depth=False, statistic="median")
    assert r["table"].index[0] == "Cell Cycle"
    with pytest.raises(ValueError, match="needs n_genes"):
        cv.occupant_contrast(*data[:4], "tumour", "normal")
    with pytest.raises(ValueError, match="statistic"):
        run(data, statistic="mode")


def test_occupants_without_a_gene_count_are_dropped_from_the_depth_balance():
    scores, near, patient, arm, genes = occupants(seed=6)
    genes = genes.copy()
    genes.iloc[::17] = np.nan
    r = cv.occupant_contrast(scores, near, patient, arm, "tumour", "normal", n_genes=genes)
    assert r["n_patients"] == 8 and np.isfinite(r["table"]["delta"]).all()


def test_no_eligible_patient_gives_an_empty_result_not_an_error():
    r = run(occupants(n_occ=4), balance_depth=False)
    assert r["n_patients"] == 0 and r["table"].empty and np.isnan(r["rho_core_signature"])


def test_the_assay_output_plugs_into_the_contrast_and_is_exported():
    """assign_cells returns the module scores and detected genes this function needs."""
    from test_pipeline_qc_aliases import make_adata
    ad = make_adata(["cytotrophoblasts", "hepatocytes"], n_per=24, seed=4)
    ad.obs["patient"] = ["P0"] * 24 + ["P1"] * 24
    ad.obs["arm"] = (["tumour"] * 12 + ["normal"] * 12) * 2
    res = cv.assign_cells(ad)
    out = pct.occupant_contrast(res["scores"], res["assignment"], ad.obs["patient"], ad.obs["arm"],
                                "tumour", "normal", region="anchor", min_cells=3, n_genes=res["n_genes"])
    assert out["table"].shape[0] == len(MODULES)


def test_benjamini_hochberg():
    p = np.array([0.01, 0.04, 0.03, 0.2, np.nan])
    q = _bh(p)
    np.testing.assert_allclose(q[:4], [0.04, 0.0533333, 0.0533333, 0.2], rtol=1e-5)
    assert np.isnan(q[4])
