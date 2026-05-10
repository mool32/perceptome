"""Validity null panels.

random_200_panel(adata, seed=42) — locked random gene panel from the dataset's own
                                  var_names; same seed + same gene order = same panel.
                                  Used to detect baseline-shift artifacts (random
                                  modules should show |log2FC| < 0.10 between any
                                  two reasonable conditions).

housekeeping_panel() — Eisenberg & Levanon 2013 minimal stable HK set.
                       Used to flag technical normalization issues.

cell_cycle_panel() — Tirosh et al. 2016 G1/S + G2/M markers.
                     Used as positive control: a real perturbation that engages
                     the cell should show |log2FC| > 0.30 on at least one cycle gene.
"""

import json
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np

_DATA = Path(__file__).parent / "data" / "validity_panels.json"


@lru_cache(maxsize=1)
def _load_panels():
    with open(_DATA) as f:
        return json.load(f)


def random_200_panel(adata, n=200, seed=42):
    """Locked random gene panel drawn from adata.var_names.

    The panel is reproducible across runs given the same seed + same adata
    var ordering. Returns list of gene symbols (length ≤ n).
    """
    rng = np.random.default_rng(seed)
    var_names = list(adata.var_names)
    if len(var_names) <= n:
        return var_names[:]
    idx = rng.choice(len(var_names), size=n, replace=False)
    return [var_names[i] for i in sorted(idx)]


def housekeeping_panel():
    return list(_load_panels()["housekeeping"]["genes"])


def cell_cycle_panel():
    return list(_load_panels()["cell_cycle"]["genes"])


def log2fc_perturbation_vs_control(
    adata,
    gene_list,
    condition_col,
    perturbation_value,
    control_value,
    pseudocount=1.0,
):
    """Mean log2 fold-change of a gene panel between two conditions.

    log2FC = log2((mean(perturbation) + pc) / (mean(control) + pc))

    Parameters
    ----------
    adata : AnnData
        Cells × genes (raw or log-normalized; method is consistent if applied to
        all conditions).
    gene_list : iterable of str
        Genes to score. Missing genes silently skipped.
    condition_col : str
    perturbation_value, control_value : str
    pseudocount : float
        Added to each side before log2 to avoid log(0).
    """
    avail = [g for g in gene_list if g in adata.var_names]
    if not avail:
        return float("nan")

    obs = adata.obs[condition_col]
    mp = obs == perturbation_value
    mc = obs == control_value
    if mp.sum() == 0 or mc.sum() == 0:
        return float("nan")

    X_p = adata[mp, avail].X
    X_c = adata[mc, avail].X
    if hasattr(X_p, "toarray"):
        X_p = X_p.toarray()
        X_c = X_c.toarray()
    mean_p = float(np.mean(X_p))
    mean_c = float(np.mean(X_c))
    return float(np.log2((mean_p + pseudocount) / (mean_c + pseudocount)))
