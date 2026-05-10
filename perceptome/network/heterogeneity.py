"""Per-module distributional analysis (Wasserstein distance)."""

import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance


def module_heterogeneity(
    scores_condition1,
    scores_condition2,
    modules=None,
):
    """Compute per-module distributional change (Wasserstein distance).

    Captures changes that mean-shift alone misses: broadening,
    bimodality, subpopulation shifts.

    Parameters
    ----------
    scores_condition1 : DataFrame
        Per-cell scores, condition 1.
    scores_condition2 : DataFrame
        Per-cell scores, condition 2.
    modules : list, optional
        Which modules to test. Default: all shared columns.

    Returns
    -------
    dict with keys:
        'wasserstein_per_module': Series (W1 distance per module)
        'mean_shift_per_module': Series (|mean1 - mean2|)
        'excess_distributional': Series (W1 - |delta_mean|, >0 = hidden heterogeneity)
        'top_heterogeneous': DataFrame (modules where excess > 0)
    """
    if modules is None:
        modules = sorted(set(scores_condition1.columns) & set(scores_condition2.columns))

    w1_vals = {}
    mean_shift = {}

    for mod in modules:
        v1 = scores_condition1[mod].dropna().values.astype(float)
        v2 = scores_condition2[mod].dropna().values.astype(float)

        if len(v1) < 10 or len(v2) < 10:
            continue

        w1 = wasserstein_distance(v1, v2)
        delta_mean = abs(np.mean(v1) - np.mean(v2))

        w1_vals[mod] = w1
        mean_shift[mod] = delta_mean

    w1_series = pd.Series(w1_vals, name='wasserstein')
    ms_series = pd.Series(mean_shift, name='mean_shift')
    excess = w1_series - ms_series
    excess.name = 'excess_distributional'

    # Modules with hidden heterogeneity
    top_het = pd.DataFrame({
        'wasserstein': w1_series,
        'mean_shift': ms_series,
        'excess': excess,
    }).sort_values('excess', ascending=False)
    top_het = top_het[top_het['excess'] > 0]

    return {
        'wasserstein_per_module': w1_series,
        'mean_shift_per_module': ms_series,
        'excess_distributional': excess,
        'top_heterogeneous': top_het,
    }
