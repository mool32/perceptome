"""Within-patient occupancy test and cross-cohort recurrence ranking.

The within-patient test is the paper's falsifiable general assay (Methods): for a cohort with
paired tumour and adjacent-normal epithelium of the same patients, the core fraction is computed
per patient for the malignant and matched-normal arms (>= 20 cells each) and compared by a
paired Wilcoxon test and a sign test. The unit of replication is the patient, not the cell.
Applied to every tissue with paired single-cell data it predicts increased occupancy with no
condition on the tissue's baseline (prostate, whose two cohorts disagree about the baseline, is
the open negative).
"""

import numpy as np
import pandas as pd
from scipy import stats

from .assign import in_region, _nearest
from .reference import load_convergence_reference, quarantined_cell_types, region_members


def _binom_p(k, n, alternative):
    if n < 1:
        return float("nan")
    return float(stats.binomtest(int(k), int(n), 0.5, alternative=alternative).pvalue)


def _wilcoxon_p(a, b, alternative):
    if len(a) < 1:
        return float("nan")
    try:
        return float(stats.wilcoxon(a, b, alternative=alternative).pvalue)
    except ValueError:        # every paired difference is zero
        return 1.0


def _matched_fraction(hit_m, genes_m, genes_n, n_draws, rng):
    """Core fraction of malignant cells resampled to the normal arm's detected-gene distribution.

    Each normal cell is paired with a malignant cell of nearest gene count (random among equally
    near ones, with replacement); the malignant fraction is averaged over ``n_draws`` draws.
    """
    candidates = []
    for g in genes_n:
        d = np.abs(genes_m - g)
        candidates.append(np.flatnonzero(d == d.min()))
    draws = [np.mean([hit_m[rng.choice(c)] for c in candidates]) for _ in range(n_draws)]
    return float(np.mean(draws))


def within_patient_occupancy(assignment, patient, arm, malignant, normal, region="core",
                             min_cells=20, n_genes=None, match_complexity=False,
                             n_draws=20, seed=0, ranking=None):
    """Paired within-patient test of region occupancy, malignant arm vs matched normal arm.

    Parameters
    ----------
    assignment : DataFrame | Series
        Output of ``assign_to_reference`` / ``assign_cells(...)['assignment']`` (or the Series
        of nearest cell types). Indexed by cell.
    patient, arm : Series
        Patient identifier and arm label per cell (aligned to ``assignment`` by index).
    malignant, normal : hashable
        Values of ``arm`` that mark the malignant and the matched-normal arm.
    region : str | int | sequence
        Region whose occupancy is compared (default: the two-cell core). See ``region_members``.
    min_cells : int
        Minimum cells per arm for a patient to enter the test (paper: 20).
    n_genes : Series | None
        Detected genes per cell. With ``match_complexity=True`` the malignant arm of each patient
        is resampled to the normal arm's gene-count distribution (``n_draws`` draws, nearest gene
        count): the control that retains ~80% of the colorectal effect in the paper (Note S-F).

    Returns
    -------
    dict
        per_patient   DataFrame: malignant / normal fraction, difference, cell counts
        n_patients, n_deeper, n_lower, n_tied
        median_malignant, median_normal, median_delta, mean_delta   (fractions, 0-1)
        wilcoxon_two_sided, wilcoxon_greater   paired Wilcoxon signed-rank P on the fractions
        sign_test_one_sided, sign_test_two_sided   binomial sign test on the sign of the difference
        matched_delta   mean complexity-matched difference (only with match_complexity)
        region, members

    Notes
    -----
    The paper reports the paired Wilcoxon two-sided and the sign test one-sided (e.g. 15 of 22
    patients: sign-test P = 0.07; 9 of 10: P = 0.011). Both sidednesses are returned.
    """
    nearest = _nearest(assignment)
    hit = in_region(assignment, region, ranking=ranking)
    df = pd.DataFrame({"hit": hit.astype(float),
                       "patient": pd.Series(patient).reindex(nearest.index),
                       "arm": pd.Series(arm).reindex(nearest.index)})
    if match_complexity:
        if n_genes is None:
            raise ValueError("match_complexity=True needs n_genes (detected genes per cell).")
        df["genes"] = pd.Series(n_genes).reindex(nearest.index).astype(float)
    df = df.dropna(subset=["patient", "arm"])

    rng = np.random.default_rng(seed)
    rows = []
    for pid, g in df.groupby("patient", sort=True):
        m, n = g[g["arm"] == malignant], g[g["arm"] == normal]
        if len(m) < min_cells or len(n) < min_cells:
            continue
        row = {"patient": pid, "malignant": m["hit"].mean(), "normal": n["hit"].mean(),
               "n_malignant": len(m), "n_normal": len(n)}
        if match_complexity:
            row["malignant_matched"] = _matched_fraction(
                m["hit"].values, m["genes"].values, n["genes"].values, n_draws, rng)
        rows.append(row)
    per = pd.DataFrame(rows)
    members = region_members(region, ranking=ranking)

    if per.empty:
        return {"per_patient": per, "n_patients": 0, "n_deeper": 0, "n_lower": 0, "n_tied": 0,
                "median_malignant": float("nan"), "median_normal": float("nan"),
                "median_delta": float("nan"), "mean_delta": float("nan"),
                "wilcoxon_two_sided": float("nan"), "wilcoxon_greater": float("nan"),
                "sign_test_one_sided": float("nan"), "sign_test_two_sided": float("nan"),
                "region": region, "members": members}

    per = per.set_index("patient")
    per["delta"] = per["malignant"] - per["normal"]
    n_deeper, n_lower = int((per["delta"] > 0).sum()), int((per["delta"] < 0).sum())
    out = {
        "per_patient": per,
        "n_patients": int(len(per)), "n_deeper": n_deeper, "n_lower": n_lower,
        "n_tied": int(len(per) - n_deeper - n_lower),
        "median_malignant": float(per["malignant"].median()), "median_normal": float(per["normal"].median()),
        "median_delta": float(per["delta"].median()), "mean_delta": float(per["delta"].mean()),
        "wilcoxon_two_sided": _wilcoxon_p(per["malignant"], per["normal"], "two-sided"),
        "wilcoxon_greater": _wilcoxon_p(per["malignant"], per["normal"], "greater"),
        "sign_test_one_sided": _binom_p(n_deeper, n_deeper + n_lower, "greater"),
        "sign_test_two_sided": _binom_p(n_deeper, n_deeper + n_lower, "two-sided"),
        "region": region, "members": members,
    }
    if match_complexity:
        per["delta_matched"] = per["malignant_matched"] - per["normal"]
        out["matched_delta"] = float(per["delta_matched"].mean())
    return out


def recurrence_ranking(assignments, top_n=None, quarantine=True):
    """Cross-cancer recurrence, depth and occupancy of each normal cell type (the paper's Table S2).

    Parameters
    ----------
    assignments : dict
        cohort name -> assignment DataFrame / Series of nearest cell types (malignant cells).
    top_n : int | None
        A cell type "recurs" in a cohort when it is among that cohort's ``top_n`` most frequent
        nearest neighbours (paper: 5). Ties are broken by name, deterministically.

    Returns
    -------
    DataFrame
        index = reference cell type; ``recurrence`` (number of cohorts), ``median_depth_pct`` and
        ``mean_depth_pct`` (per-cohort fraction of cells reaching it), ``occupancy_pct`` (pooled
        over all cells). Sorted by recurrence, then median and mean depth.
    """
    ref = load_convergence_reference()
    if top_n is None:
        top_n = ref["definitions"]["qc"]["recurrence_top_n"]
    drop = set(quarantined_cell_types()) if quarantine else set()
    names = [c for c in ref["tables"]["S2"].index if c not in drop]

    fractions, counts = {}, pd.Series(0, index=names, dtype=float)
    recur = pd.Series(0, index=names, dtype=int)
    for cohort, a in assignments.items():
        near = _nearest(a).dropna()
        if len(near) == 0:
            continue
        c = near.value_counts().reindex(names).fillna(0)
        fractions[cohort] = c / len(near)
        counts += c
        order = sorted(names, key=lambda n: (-c[n], n))
        recur[order[:top_n]] += 1
    frac = pd.DataFrame(fractions) if fractions else pd.DataFrame(index=names)
    total = counts.sum()
    out = pd.DataFrame({
        "recurrence": recur,
        "median_depth_pct": 100 * frac.median(axis=1) if len(frac.columns) else 0.0,
        "mean_depth_pct": 100 * frac.mean(axis=1) if len(frac.columns) else 0.0,
        "occupancy_pct": 100 * counts / total if total else 0.0,
    })
    return out.sort_values(["recurrence", "median_depth_pct", "mean_depth_pct"], ascending=False, kind="mergesort")
