"""What malignancy adds inside the shared state: malignant versus normal occupants (paper Fig. 6, Fig. S18).

Increased occupancy is only half of what malignancy does. The sharper question is how a malignant
cell already in the convergent state differs from a normal cell of the same patient that is also in
it — same state, same tissue, same patient, so lineage and batch cancel. The paper compares the two
occupant populations module by module within each of thirty colorectal patients (>= 10 core-assigned
cells per arm), balancing detection depth so the comparison is not a library-size effect (median
2,640 versus 2,496 detected genes, from a raw 4,351 versus 2,395); sixteen modules separate the
occupants reproducibly (stress sensors and anabolic programs up in the malignant occupant; PPARalpha,
MR and PI3K/PTEN down), and the differences correlate only weakly with the axis that defines the
core (Spearman rho = +0.36 over the 44 modules): the difference lies largely off that axis. The paper
presents this for colorectal only ("awaits replication").

The manuscript gives the unit (the patient), the cell floor, the paired Wilcoxon test and the depth
balance, but not the rule that selects the sixteen modules; ``occupant_contrast`` therefore returns the
per-module statistics and leaves the threshold to the caller.
"""

import numpy as np
import pandas as pd
from scipy import stats

from .assign import in_region, _nearest
from .occupancy import _wilcoxon_p, _nearest_candidates, _draw_matches
from .signature import core_signature


def _bh(p):
    """Benjamini-Hochberg adjusted P values."""
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order] * len(q) / (np.arange(len(q)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(len(q))
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def occupant_contrast(scores, assignment, patient, arm, malignant, normal, region="core",
                      min_cells=10, n_genes=None, balance_depth=True, n_draws=20, seed=0,
                      statistic="mean", calibration="pseudobulk", ranking=None):
    """Per-module difference between malignant and matched-normal occupants of a region.

    Parameters
    ----------
    scores : DataFrame
        Module scores per cell (``assign_cells(...)['scores']``).
    assignment : DataFrame | Series
        Nearest normal cell type per cell (``assign_cells(...)['assignment']``); the *occupants* are
        the cells assigned to ``region``.
    patient, arm : Series
        Patient and arm label per cell; ``malignant`` and ``normal`` are the values of ``arm``.
    region : str | int | sequence
        The state whose occupants are compared (default: the two-cell core).
    min_cells : int
        Minimum occupants per arm for a patient to enter (paper: 10).
    n_genes : Series | None
        Detected genes per cell (``assign_cells(...)['n_genes']``). With ``balance_depth`` the
        malignant occupants of each patient are resampled to the normal occupants' gene-count
        distribution (nearest gene count, ``n_draws`` draws, as in Note S-F), so that a depth
        difference between the arms cannot drive the contrast.
    statistic : 'mean' | 'median'
        Per-patient summary of a module over an arm's occupants.
    calibration : str | Calibration
        Calibration that defines the core signature for the comparison with the defining axis.

    Returns
    -------
    dict
        table        per module: delta_raw (mean over patients of the uncontrolled malignant-minus-
                     normal difference), delta (the depth-balanced difference the test uses; equal to
                     delta_raw without balancing), n_higher (patients with delta > 0), wilcoxon_p
                     (paired, two-sided, over patients), q_bh (Benjamini-Hochberg over the modules —
                     an addition of this package, the paper states no multiplicity rule). Sorted by P.
        per_patient  patients x modules matrix of the balanced differences
        n_cells      occupants per patient and arm
        n_patients, patients
        rho_core_signature, p_core_signature   Spearman correlation over the modules between ``delta``
                     and the standardized core signature (paper: +0.36, colorectal)
        depth        with ``n_genes``: median detected genes of malignant and normal occupants, raw and
                     after balancing (paper: 4,351 vs 2,395 raw, 2,640 vs 2,496 balanced)
    """
    if statistic not in ("mean", "median"):
        raise ValueError("statistic must be 'mean' or 'median'.")
    if balance_depth and n_genes is None:
        raise ValueError("balance_depth=True needs n_genes (detected genes per cell); pass "
                         "balance_depth=False to compare without the depth control.")
    near = _nearest(assignment)
    occ = in_region(assignment, region, ranking=ranking).reindex(near.index).fillna(False).astype(bool)
    S = pd.DataFrame(scores).reindex(near.index)
    info = pd.DataFrame({"patient": pd.Series(patient).reindex(near.index),
                         "arm": pd.Series(arm).reindex(near.index), "occ": occ})
    if n_genes is not None:
        info["genes"] = pd.Series(n_genes).reindex(near.index).astype(float)
    info = info[info["occ"] & info["patient"].notna() & info["arm"].notna()]
    if balance_depth:
        info = info[info["genes"].notna()]
    summarize = (lambda a: a.mean(axis=0)) if statistic == "mean" else (lambda a: np.median(a, axis=0))

    rng = np.random.default_rng(seed)
    raw_rows, ctl_rows, counts, patients = [], [], [], []
    depth_raw_m, depth_raw_n, depth_bal_m = [], [], []
    for pid, g in info.groupby("patient", sort=True):
        m, n = g[g["arm"] == malignant], g[g["arm"] == normal]
        if len(m) < min_cells or len(n) < min_cells:
            continue
        Sm, Sn = S.loc[m.index].values, S.loc[n.index].values
        normal_stat = summarize(Sn)
        raw = summarize(Sm) - normal_stat
        if balance_depth:
            cand = _nearest_candidates(m["genes"].values, n["genes"].values)
            picks = _draw_matches(cand, n_draws, rng)
            ctl = np.mean([summarize(Sm[p]) for p in picks], axis=0) - normal_stat
            depth_bal_m.append(np.concatenate([m["genes"].values[p] for p in picks]))
        else:
            ctl = raw
        patients.append(pid)
        raw_rows.append(raw)
        ctl_rows.append(ctl)
        counts.append({"patient": pid, "n_malignant": len(m), "n_normal": len(n)})
        if n_genes is not None:
            depth_raw_m.append(m["genes"].values)
            depth_raw_n.append(n["genes"].values)

    modules = list(S.columns)
    if not patients:
        empty = pd.DataFrame(columns=["delta_raw", "delta", "n_higher", "wilcoxon_p", "q_bh"])
        return {"table": empty, "per_patient": pd.DataFrame(columns=modules), "n_cells": pd.DataFrame(),
                "n_patients": 0, "patients": [], "rho_core_signature": float("nan"),
                "p_core_signature": float("nan"), "depth": None}

    D = pd.DataFrame(np.array(ctl_rows), index=patients, columns=modules)
    Draw = pd.DataFrame(np.array(raw_rows), index=patients, columns=modules)
    pvals = [_wilcoxon_p(D[mod].values, np.zeros(len(D)), "two-sided") for mod in modules]
    table = pd.DataFrame({
        "delta_raw": Draw.mean(axis=0), "delta": D.mean(axis=0),
        "n_higher": (D > 0).sum(axis=0).astype(int), "wilcoxon_p": pvals, "q_bh": _bh(pvals),
    }, index=modules)
    table = table.sort_values(["wilcoxon_p", "delta"], ascending=[True, False], kind="mergesort")

    sig = core_signature(calibration, region).reindex(modules)
    both = pd.concat([table["delta"], sig], axis=1, keys=["d", "s"]).dropna()
    if len(both) > 2 and both["d"].nunique() > 1:
        rho, rp = stats.spearmanr(both["d"], both["s"])
    else:
        rho, rp = float("nan"), float("nan")

    depth = None
    if n_genes is not None:
        raw_m, raw_n = np.concatenate(depth_raw_m), np.concatenate(depth_raw_n)
        depth = {"malignant_raw_median": float(np.median(raw_m)), "normal_raw_median": float(np.median(raw_n))}
        if balance_depth:
            depth["malignant_balanced_median"] = float(np.median(np.concatenate(depth_bal_m)))
            depth["normal_balanced_median"] = depth["normal_raw_median"]
    return {"table": table, "per_patient": D, "n_cells": pd.DataFrame(counts).set_index("patient"),
            "n_patients": int(len(patients)), "patients": patients,
            "rho_core_signature": float(rho), "p_core_signature": float(rp), "depth": depth}
