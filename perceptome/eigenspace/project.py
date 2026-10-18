"""Project module scores into the fixed 9-PC reference eigenspace.

The paper's coordinate system is *standardized*: module scores are z-scored
against the reference calibration and then projected through the fixed loading
matrix (Methods; Note S1). ``project`` follows that by default. Earlier tool
versions (<= 0.2.2) multiplied raw scores by the loadings, which leaves every
cell sharing one large constant offset (mean pairwise cosine between the 154
normal cell types 0.97, instead of ~0 in the standardized space); that behaviour
is still available as ``standardize=False`` and is not the paper's geometry.
"""

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import Calibration, load_calibration, warn_if_scale_mismatch

_DEFAULT_REF = Path(__file__).parent / "data" / "reference_v03.json"


@lru_cache(maxsize=2)
def _load_reference(path=None):
    p = Path(path) if path else _DEFAULT_REF
    if not p.exists():
        raise FileNotFoundError(
            f"Eigenspace reference not found at {p}. "
            "Run scripts/03_build_eigenspace.py to generate."
        )
    with open(p) as f:
        return json.load(f)


def _resolve_reference(reference):
    return _load_reference(reference) if (reference is None or isinstance(reference, str)) else reference


def _resolve_calibration(ref, calibration):
    """Calibration to standardize with: explicit one, a named one, or the reference's own stats."""
    if calibration is None:
        return Calibration.from_reference(ref)
    if isinstance(calibration, Calibration):
        return calibration
    if isinstance(calibration, str):
        return load_calibration(calibration, reference=ref)
    raise TypeError(f"calibration must be None, a name or a Calibration, got {type(calibration).__name__}")


def project(scores, reference=None, n_pcs=None, confidence_flags=True,
            standardize=True, calibration=None, missing_modules=None, warn_scale=True):
    """Project module scores into the 9-PC reference eigenspace.

    Parameters
    ----------
    scores : DataFrame
        cells/clusters x modules (``mean_raw`` scores). Modules that are not columns
        are treated as unmeasured (placed at the reference mean).
    reference : str | dict | None
        Eigenspace reference JSON path or pre-loaded dict (loadings).
    n_pcs : int | None
        Number of PCs to return (default: all PCs of the reference).
    confidence_flags : bool
        Attach per-PC bootstrap-stability label.
    standardize : bool
        True (default): z-score each module against the calibration, then project —
        the paper's coordinate system. False: legacy raw projection (``scores @ loadings``).
    calibration : None | 'pseudobulk' | 'single_cell_scaled' | Calibration
        Reference scaling used for z-scoring. Default: the pseudobulk statistics stored in
        the eigenspace reference. Per-cell data need ``single_cell_scaled`` (see
        :mod:`perceptome.eigenspace.calibration`).
    missing_modules : iterable of str | None
        Modules to treat as unmeasured even though a column exists (e.g. modules whose genes
        are absent from the dataset: ``coverage.index[coverage.n_genes_found < 2]``).
    warn_scale : bool
        Warn when the scores look like per-cell values under a pseudobulk calibration.

    Returns
    -------
    dict
        coordinates           DataFrame (samples x n_pcs)
        stable_coordinates    DataFrame (samples x top stable PCs only)
        module_contributions  DataFrame (modules x n_pcs) — the loadings
        confidence            dict {PC: 'stable'|'probable'|'exploratory'}
        calibration           name of the calibration used, or None when standardize=False
    """
    ref = _resolve_reference(reference)
    module_order = ref["module_order"]
    n_pcs = ref["n_pcs"] if n_pcs is None else min(int(n_pcs), ref["n_pcs"])
    pc_names = [f"PC{i + 1}" for i in range(n_pcs)]

    loadings = np.zeros((len(module_order), n_pcs))
    for j in range(n_pcs):
        pc_loads = ref["loadings"][f"PC{j + 1}"]
        for i, mod in enumerate(module_order):
            loadings[i, j] = pc_loads.get(mod, 0.0)

    if standardize:
        cal = _resolve_calibration(ref, calibration)
        Z = cal.standardize(scores)[module_order]
        if missing_modules:
            Z.loc[:, [m for m in missing_modules if m in Z.columns]] = 0.0
        if warn_scale:
            warn_if_scale_mismatch(Z, cal)
        coords = Z.values @ loadings
        cal_name = cal.name
    else:
        aligned = np.zeros((len(scores), len(module_order)))
        for i, mod in enumerate(module_order):
            if mod in scores.columns:
                aligned[:, i] = scores[mod].values
        coords = aligned @ loadings
        cal_name = None

    coords_df = pd.DataFrame(coords, index=scores.index, columns=pc_names)
    contrib_df = pd.DataFrame(loadings, index=module_order, columns=pc_names)

    confidence = {}
    if confidence_flags:
        stability = ref.get("bootstrap_stability", [])
        for i, pc in enumerate(pc_names):
            if pc in ref.get("pc_confidence", {}):
                confidence[pc] = ref["pc_confidence"][pc]
            elif i < len(stability):
                s = stability[i]
                confidence[pc] = "stable" if s > 0.6 else "probable" if s > 0.4 else "exploratory"
            else:
                confidence[pc] = "exploratory"

    stable_pcs = [pc for pc, conf in confidence.items() if conf == "stable"] or pc_names[:3]
    stable_df = coords_df[stable_pcs].copy()

    return {
        "coordinates": coords_df,
        "stable_coordinates": stable_df,
        "module_contributions": contrib_df,
        "confidence": confidence,
        "calibration": cal_name,
    }


def project_shift(delta, reference=None, n_pcs=None, calibration=None):
    """Project a module-score *difference* (tumour - normal, drug - control, ...).

    A shift is carried into the standardized coordinate system by dividing each module
    by its reference SD and projecting through the loadings (no centering), so that
    ``project(a)['coordinates'] - project(b)['coordinates'] == project_shift(a - b)``.

    Parameters
    ----------
    delta : Series (modules) | DataFrame (shifts x modules)
    """
    ref = _resolve_reference(reference)
    cal = _resolve_calibration(ref, calibration)
    return cal.project_shift(delta, n_pcs=n_pcs)
