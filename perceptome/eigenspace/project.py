"""Project module scores into the canonical 12-PC eigenspace."""

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

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


def project(scores, reference=None, n_pcs=12, confidence_flags=True):
    """Project module scores into the canonical 12-PC eigenspace.

    Parameters
    ----------
    scores : DataFrame
        cells/clusters × modules. Module name strings used; missing modules → 0.
    reference : str | dict | None
        Reference eigenspace JSON path or pre-loaded dict.
    n_pcs : int
        How many PCs to project onto (max = reference n_pcs).
    confidence_flags : bool
        Attach per-PC bootstrap-stability label.

    Returns
    -------
    dict
        coordinates           DataFrame (samples × n_pcs)
        stable_coordinates    DataFrame (samples × top stable PCs only)
        module_contributions  DataFrame (modules × n_pcs)
        confidence            dict {PC: 'stable'|'probable'|'exploratory'}
    """
    ref = _load_reference(reference) if (reference is None or isinstance(reference, str)) else reference

    module_order = ref["module_order"]
    n_pcs = min(n_pcs, ref["n_pcs"])

    loadings = np.zeros((len(module_order), n_pcs))
    for j in range(n_pcs):
        pc_loads = ref["loadings"][f"PC{j+1}"]
        for i, mod in enumerate(module_order):
            loadings[i, j] = pc_loads.get(mod, 0.0)

    aligned = np.zeros((len(scores), len(module_order)))
    for i, mod in enumerate(module_order):
        if mod in scores.columns:
            aligned[:, i] = scores[mod].values

    coords = aligned @ loadings
    pc_names = [f"PC{i+1}" for i in range(n_pcs)]
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
    }
