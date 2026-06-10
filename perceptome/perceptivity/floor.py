"""Capacity-floor predictor (Factor 1) + operation-intensity combination (Factor 2).

Closed across paper4.5 + paper4.6 + paper4.7 + paper4.8 (2026-05-10):
  - A_baseline > 4.5  ⇒  saturated, no upward ramp possible
  - A_baseline < 2.5  ⇒  capacious, ramp possible (magnitude operation-determined)
  - 2.5 ≤ A ≤ 4.5    ⇒  intermediate
  - Downward suppression (specific signaling, e.g. atRA → UPR-ATF6) is allowed —
    the predictor is upward-asymmetric, not absolute.

Two-factor framework (v0.4 Consolidation):
    engagement(cell, op, module) = capacity(cell, module) × operation_intensity(op, module)
  Factor 1 (capacity) — HPA-derivable, always returned.
  Factor 2 (operation_intensity) — calibrated from the substrate-series corpus
    (operation_intensity_v1), returned ONLY when `operation` is supplied. Ordinal
    bands, NOT a point estimate (data-limited: ~10 paradigms, strict 3a FAILed).
"""

from typing import Iterable, Sequence

import numpy as np
import pandas as pd


SATURATED_BLOCKED_UP = "saturated_blocked_up"
CAPACIOUS = "capacious"
INTERMEDIATE = "intermediate"
NO_DATA = "no_data"

# operation taxonomy calibrated in operation_intensity_v1 (substrate-series 4.3-4.8 + paper4)
OPERATIONS = (
    "immune_activation",
    "hypertrophy",
    "terminal_differentiation",
    "differentiation",
    "regeneration",
    "retinoid_perturbation",
    "memory_consolidation",
)


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


def _combine_factors(floor_label, intens):
    """Combine Factor 1 (capacity floor) × Factor 2 (operation intensity).

    Returns (predicted_direction, predicted_magnitude). `intens` is an
    operation_intensity_v1 record, or None when the operation has no
    calibration for this module.
    """
    if intens is None:
        return "unknown", "unknown"  # no Factor-2 calibration for this (op, module)
    band = intens["intensity_band"]
    idir = intens["direction"]
    # robust = |effect| ≥ 0.30 (medium+) — matches the substrate papers' own
    # engagement threshold; below it, effects sit within several studies' null
    # floor (random-200 ran up to 0.6 d in paper4.5), so direction is not trusted.
    robust = band in ("medium", "high", "extreme")
    if band == "none" or idir == "flat":
        return "flat", "none"                     # operation does not engage this module
    if idir == "suppressive":
        # upward-asymmetric: active suppression allowed regardless of capacity,
        # but a confident 'down' needs a robust downward push (e.g. atRA→UPR-ATF6).
        return ("down", band) if robust else ("down_weak", band)
    if floor_label == SATURATED_BLOCKED_UP:
        return "up_blocked", "blocked"            # Factor-1 ceiling wins
    if floor_label == CAPACIOUS:
        return "up", band                         # operation can ramp it
    if floor_label == INTERMEDIATE:
        return "up_uncertain", f"attenuated:{band}"
    return "unknown", band                        # capacity unknown


def predict_engagement(
    starting_cell_type,
    operation_modules=None,
    operation=None,
    hpa_reference=None,
    intensity_reference=None,
    hi=4.5,
    lo=2.5,
):
    """Predict architecture engagement for a cell × (operation) × module.

    Two modes:

    * **Factor 1 only** (``operation=None``) — capacity-floor classification per
      module from HPA. ``predicted_magnitude`` is ``'unknown'`` (no operation
      context). Backward-compatible with v0.2.
    * **Two-factor** (``operation`` given) — combines the cell's capacity floor
      with the calibrated operation_intensity (Factor 2) to yield an ordinal
      ``predicted_magnitude`` and a direction (incl. ``down`` for operations that
      actively suppress a module, e.g. retinoid → UPR-ATF6).

    Parameters
    ----------
    starting_cell_type : str
        Cell type name (matched against HPA reference, case-insensitive substring).
    operation_modules : str | Iterable[str] | None
        Module(s) to predict. If None and ``operation`` is given, defaults to all
        modules calibrated for that operation. Required when ``operation`` is None.
    operation : str | None
        One of ``OPERATIONS``. Enables Factor 2.
    hpa_reference, intensity_reference : dict | None
        Loaded automatically if None.
    hi, lo : float
        Capacity-floor thresholds.

    Returns
    -------
    DataFrame indexed by module. Factor-1 columns always present:
        R_baseline, A_baseline, C, headroom, capacity_floor, predicted_direction,
        predicted_magnitude.
    Factor-2 mode adds:
        operation, operation_intensity, intensity_direction, intensity_metric,
        n_anchors, confidence.
    """
    from .reference import load_hpa_perceptivity, load_operation_intensity

    if hpa_reference is None:
        hpa_reference = load_hpa_perceptivity()
    R = hpa_reference["R"]
    A = hpa_reference["A"]
    C = hpa_reference["C"]
    H = hpa_reference["headroom"]

    itab = None
    if operation is not None:
        if operation not in OPERATIONS:
            raise ValueError(f"operation must be one of {OPERATIONS}, got {operation!r}")
        if intensity_reference is None:
            intensity_reference = load_operation_intensity()
        itab = intensity_reference["table"]
        if operation_modules is None:
            operation_modules = sorted({m for (op, m) in itab if op == operation})
    elif operation_modules is None:
        raise ValueError(
            "provide operation_modules (Factor-1 only) or operation (Factor-2 two-factor mode)"
        )

    if isinstance(operation_modules, str):
        operation_modules = [operation_modules]
    operation_modules = list(operation_modules)

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
            row = {
                "module": mod, "R_baseline": np.nan, "A_baseline": np.nan,
                "C": np.nan, "headroom": np.nan,
                "capacity_floor": NO_DATA,
                "predicted_direction": "unknown_module_not_in_catalog",
                "predicted_magnitude": "unknown",
            }
            if operation is not None:
                row.update(operation=operation, operation_intensity="no_data",
                           intensity_direction="no_data", intensity_metric=None,
                           n_anchors=0, confidence="none")
            rows.append(row)
            continue

        a = float(A.loc[starting_cell_type, mod])
        floor = capacity_floor(a, hi=hi, lo=lo)

        if operation is None:
            # Factor-1-only (backward compatible)
            if floor == SATURATED_BLOCKED_UP:
                direction = "up_blocked_down_possible"
            elif floor == CAPACIOUS:
                direction = "up_possible"
            elif floor == INTERMEDIATE:
                direction = "intermediate_uncertain"
            else:
                direction = "unknown"
            magnitude = "unknown"
            extra = {}
        else:
            intens = itab.get((operation, mod))
            direction, magnitude = _combine_factors(floor, intens)
            if intens is None:
                extra = dict(operation=operation, operation_intensity="no_data",
                             intensity_direction="no_data", intensity_metric=None,
                             n_anchors=0, confidence="none")
            else:
                low_conf = intens.get("n_obs", 0) < 2 or bool(intens.get("calib_saturated"))
                extra = dict(operation=operation,
                             operation_intensity=intens["intensity_band"],
                             intensity_direction=intens["direction"],
                             intensity_metric=intens.get("effect_metric"),
                             n_anchors=int(intens.get("n_obs", 0)),
                             confidence="low" if low_conf else "medium")

        row = {
            "module": mod,
            "R_baseline": float(R.loc[starting_cell_type, mod]),
            "A_baseline": a,
            "C": float(C.loc[starting_cell_type, mod]),
            "headroom": float(H.loc[starting_cell_type, mod]),
            "capacity_floor": floor,
            "predicted_direction": direction,
            "predicted_magnitude": magnitude,
        }
        row.update(extra)
        rows.append(row)

    return pd.DataFrame(rows).set_index("module")
