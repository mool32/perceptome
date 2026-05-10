"""Compute the 5-vector perceptivity metric per cell type × module.

Inputs are R and A score matrices (cells × modules) plus an optional cell-type
column and tissue/class column. Output is a per-cell-type DataFrame with the
five components plus the categorical specialization quadrant.

Per the v1.4 amendment, R/A/C/headroom are continuous (not binary). Binary
quantities (BS/GS engagement count, quadrant) use an absolute A threshold,
default 2.0 (matches v1.3/v1.4 thresholds). All thresholds are user-overridable.
"""

import numpy as np
import pandas as pd

from .floor import capacity_floor


_DEFAULT_A_THRESHOLD = 2.0  # absolute log1p(nCPM) threshold for "engaged" module
_DEFAULT_FLOOR_HI = 4.5     # A_baseline > 4.5 ⇒ saturated (per v1.4)
_DEFAULT_FLOOR_LO = 2.5     # A_baseline < 2.5 ⇒ capacious (per v1.4)


def compute_perceptivity(
    R_scores,
    A_scores,
    cell_type=None,
    cell_class=None,
    A_threshold=_DEFAULT_A_THRESHOLD,
    floor_hi=_DEFAULT_FLOOR_HI,
    floor_lo=_DEFAULT_FLOOR_LO,
    A_max=None,
):
    """Compute the 5-vector perceptivity metric per (cell type, module) and aggregate per cell type.

    Parameters
    ----------
    R_scores, A_scores : DataFrame
        cells × modules. Same index, same columns. Both are mean log expression
        (`mean_raw` from score_modules with gene_set='core' / gene_set='activity').
    cell_type : Series | None
        adata.obs cell type column (index aligned to R/A). If None, returns
        per-cell perceptivity (no aggregation).
    cell_class : Series | None
        Tissue / cell-class column for BS computation. Index aligned to R/A.
        If None, BS is set to NaN (only GS computed).
    A_threshold : float
        Absolute A above which a module is considered "engaged". Default 2.0.
    floor_hi, floor_lo : float
        Saturation / capacity thresholds for the capacity-floor classifier.
    A_max : Series | None
        A_max per module (for headroom). If None, computed from the input
        as the max across cell types — but for cross-dataset analysis you
        usually want the HPA-derived A_max from load_hpa_perceptivity().

    Returns
    -------
    dict
        per_cell_type      DataFrame (cell_types × {R, A, C, I, BS, GS, spec_quadrant, n_engaged})
        per_module         dict {component: DataFrame (cell_types × modules)}
                           keys: R, A, C, headroom, A_engaged_mask, capacity_floor
    """
    if list(R_scores.columns) != list(A_scores.columns):
        raise ValueError("R_scores and A_scores must have identical column order")
    if list(R_scores.index) != list(A_scores.index):
        raise ValueError("R_scores and A_scores must have identical row index")

    if cell_type is None:
        # Treat each row as its own "cell type"
        ct = pd.Series(R_scores.index, index=R_scores.index)
    else:
        ct = pd.Series(cell_type, index=R_scores.index)

    R_ct = R_scores.groupby(ct).mean()
    A_ct = A_scores.groupby(ct).mean()
    C_ct = R_ct - A_ct

    if A_max is None:
        A_max = A_ct.max(axis=0)
    else:
        A_max = pd.Series(A_max).reindex(R_ct.columns)
    headroom_ct = A_max - A_ct

    engaged = (A_ct > A_threshold).astype(int)
    floor_ct = R_ct.copy()
    for col in R_ct.columns:
        floor_ct[col] = [capacity_floor(a, hi=floor_hi, lo=floor_lo) for a in A_ct[col].values]

    # Per-cell-type aggregates
    R_agg = engaged.sum(axis=1)
    A_agg = engaged.sum(axis=1)  # number of engaged modules; same as R_binary in v1.4 sense
    # Note: in v1.4 R_agg used R-binary (R > 2.0). Here we keep the more directly meaningful
    # "n_engaged_activity" for both R/A counts in the aggregate; the per-module R/A are
    # still in per_module['R'] / per_module['A'] for users who need the readiness side.

    I_mean = A_ct.where(engaged == 1).mean(axis=1)
    n_engaged = engaged.sum(axis=1)

    BS, GS = _compute_specialization_breadth(engaged, cell_class)

    bs_med = BS.median(skipna=True)
    gs_med = GS.median(skipna=True)
    quadrant = pd.Series(index=engaged.index, dtype=object)
    for ct_name in engaged.index:
        bs = BS.get(ct_name, np.nan)
        gs = GS.get(ct_name, np.nan)
        quadrant[ct_name] = classify_quadrant(bs, gs, bs_med, gs_med)

    per_ct = pd.DataFrame({
        "n_engaged": n_engaged,
        "I_mean": I_mean,
        "BS": BS,
        "GS": GS,
        "spec_quadrant": quadrant,
    })

    return {
        "per_cell_type": per_ct,
        "per_module": {
            "R": R_ct,
            "A": A_ct,
            "C": C_ct,
            "headroom": headroom_ct,
            "engaged_mask": engaged,
            "capacity_floor": floor_ct,
        },
    }


def perceptivity_per_celltype(
    adata,
    cell_type_column,
    cell_class_column=None,
    score_method="mean_raw",
    catalog=None,
    A_threshold=_DEFAULT_A_THRESHOLD,
):
    """Convenience wrapper: scores R and A from adata, then computes perceptivity.

    Parameters
    ----------
    adata : AnnData
        Cells × genes (log-normalized).
    cell_type_column : str
        adata.obs column for cell-type aggregation.
    cell_class_column : str | None
        adata.obs column for tissue/class breadth (for BS).
    score_method : 'mean_raw' | 'scanpy_corrected'
        Scoring method. mean_raw is the default and matches the v1.4 metric.
    """
    from ..score import score_modules

    R = score_modules(adata, gene_set="core", method=score_method, catalog=catalog)["scores"]
    A = score_modules(adata, gene_set="activity", method=score_method, catalog=catalog)["scores"]

    ct = adata.obs[cell_type_column]
    klass = adata.obs[cell_class_column] if cell_class_column else None

    return compute_perceptivity(
        R, A, cell_type=ct, cell_class=klass, A_threshold=A_threshold,
    )


def classify_quadrant(BS, GS, BS_median, GS_median):
    """Categorical 4-quadrant typology from BS / GS relative to dataset medians.

    Q1 unique         high BS, high GS
    Q2 locally-homogeneous-globally-unique-class   low BS, high GS
    Q3 locally-divided-labor-globally-common       high BS, low GS
    Q4 generalist     low BS, low GS

    Returns string label or 'NA' if either is NaN.
    """
    if BS is None or GS is None or BS != BS or GS != GS:  # NaN check
        return "NA"
    bh = BS >= BS_median
    gh = GS >= GS_median
    if bh and gh:
        return "Q1: locally + globally unique"
    if (not bh) and gh:
        return "Q2: locally homogeneous, globally unique class"
    if bh and (not gh):
        return "Q3: locally divided-labor, globally common"
    return "Q4: generalist"


def _compute_specialization_breadth(engaged_df, cell_class):
    """BS and GS per cell type from engagement matrix.

    For each cell type C with set of engaged modules E_C:
      BS(C) = sum over m in E_C of (1 - frac_engaged_in(m, neighbors_in_class))
      GS(C) = sum over m in E_C of (1 - frac_engaged_in(m, all_other_cells))
    """
    cell_types = engaged_df.index
    BS = pd.Series(index=cell_types, dtype=float)
    GS = pd.Series(index=cell_types, dtype=float)

    if cell_class is not None:
        klass = pd.Series(cell_class).reindex(cell_types).fillna("__unknown__")
    else:
        klass = None

    for ct in cell_types:
        engaged_modules = engaged_df.columns[engaged_df.loc[ct].astype(bool)]
        if len(engaged_modules) == 0:
            BS[ct] = 0.0
            GS[ct] = 0.0
            continue

        others_global = [c for c in cell_types if c != ct]
        if others_global:
            frac_g = engaged_df.loc[others_global, engaged_modules].mean(axis=0)
            GS[ct] = float((1 - frac_g).sum())
        else:
            GS[ct] = float("nan")

        if klass is not None:
            ct_klass = klass[ct]
            neighbors = [c for c in others_global if klass.get(c) == ct_klass]
            if neighbors:
                frac_b = engaged_df.loc[neighbors, engaged_modules].mean(axis=0)
                BS[ct] = float((1 - frac_b).sum())
            else:
                BS[ct] = float("nan")
        else:
            BS[ct] = float("nan")

    return BS, GS
