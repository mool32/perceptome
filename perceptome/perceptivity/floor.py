"""Capacity-floor predictor (upward-asymmetric).

Closed across paper4.5 + paper4.6 + paper4.7 + paper4.8 (2026-05-10):
  - A_baseline > 4.5  ⇒  saturated, no upward ramp possible
  - A_baseline < 2.5  ⇒  capacious, ramp possible (magnitude operation-determined)
  - 2.5 ≤ A ≤ 4.5    ⇒  intermediate
  - Downward suppression (specific signaling, e.g. atRA → UPR-ATF6) is allowed —
    the predictor is upward-asymmetric, not absolute.

predict_engagement() applies the two-factor framework: returns the HPA-derivable
factor 1 (capacity), and a qualitative direction prediction. Magnitude requires
operation-intensity modeling outside the tool.
"""

from typing import Iterable, Sequence

import numpy as np
import pandas as pd


SATURATED_BLOCKED_UP = "saturated_blocked_up"
CAPACIOUS = "capacious"
INTERMEDIATE = "intermediate"
NO_DATA = "no_data"


def capacity_floor(A_baseline, hi=4.5, lo=2.5):
    """Classify a single A_baseline into capacity-floor regime.

    Returns one of: 'saturated_blocked_up' | 'capacious' | 'intermediate' | 'no_data'.
    """
    if A_baseline is None or (isinstance(A_baseline, float) and np.isnan(A_baseline)):
        return NO_DATA
    if A_baseline > hi:
        return SATURATED_BLOCKED_UP
    if A_baseline < lo:
        return CAPACIOUS
    return INTERMEDIATE


def predict_engagement(
    starting_cell_type,
    operation_modules,
    hpa_reference=None,
    hi=4.5,
    lo=2.5,
):
    """Predict architecture engagement for a cell × operation pair, factor 1 only.

    Per the v1.4 two-factor framework:
        engagement(cell, op, module) = capacity(cell, module) × operation_intensity(op, module)

    Tool covers Factor 1 (capacity) only — derivable from HPA. Factor 2
    (operation intensity) requires context modeling and is returned as 'unknown'.

    Parameters
    ----------
    starting_cell_type : str
        Cell type name (must match an HPA cell type or user-supplied reference).
    operation_modules : str | Iterable[str]
        Module(s) that the operation activates.
    hpa_reference : dict | None
        From load_hpa_perceptivity(). If None, loaded automatically.
    hi, lo : float
        Capacity-floor thresholds.

    Returns
    -------
    DataFrame
        One row per requested module with columns:
          R_baseline, A_baseline, C, headroom, capacity_floor,
          predicted_direction (up_blocked | up_possible | down_only | unknown),
          predicted_magnitude (always 'unknown' until factor 2 is modeled).
    """
    from .reference import load_hpa_perceptivity

    if isinstance(operation_modules, str):
        operation_modules = [operation_modules]
    operation_modules = list(operation_modules)

    if hpa_reference is None:
        hpa_reference = load_hpa_perceptivity()

    R = hpa_reference["R"]
    A = hpa_reference["A"]
    C = hpa_reference["C"]
    H = hpa_reference["headroom"]

    if starting_cell_type not in R.index:
        candidates = [c for c in R.index if starting_cell_type.lower() in c.lower()]
        if not candidates:
            raise KeyError(
                f"Cell type {starting_cell_type!r} not in HPA reference. "
                f"Try one of {sorted(R.index)[:5]}... ({len(R.index)} total)."
            )
        starting_cell_type = candidates[0]

    rows = []
    for mod in operation_modules:
        if mod not in R.columns:
            rows.append({
                "module": mod, "R_baseline": np.nan, "A_baseline": np.nan,
                "C": np.nan, "headroom": np.nan,
                "capacity_floor": NO_DATA,
                "predicted_direction": "unknown_module_not_in_catalog",
                "predicted_magnitude": "unknown",
            })
            continue
        a = float(A.loc[starting_cell_type, mod])
        floor = capacity_floor(a, hi=hi, lo=lo)
        if floor == SATURATED_BLOCKED_UP:
            direction = "up_blocked_down_possible"
        elif floor == CAPACIOUS:
            direction = "up_possible"
        elif floor == INTERMEDIATE:
            direction = "intermediate_uncertain"
        else:
            direction = "unknown"
        rows.append({
            "module": mod,
            "R_baseline": float(R.loc[starting_cell_type, mod]),
            "A_baseline": a,
            "C": float(C.loc[starting_cell_type, mod]),
            "headroom": float(H.loc[starting_cell_type, mod]),
            "capacity_floor": floor,
            "predicted_direction": direction,
            "predicted_magnitude": "unknown",  # factor 2 not modeled
        })

    return pd.DataFrame(rows).set_index("module")
