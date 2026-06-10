"""Bundled HPA perceptivity reference (154 cell types × 44 modules).

Precomputed from HPA single-cell-type RNA expression (rna_single_cell_type.tsv,
~3M rows) using mean log1p(nCPM) of core_genes (R) and activity_genes (A) per
(cell type, module). C = R − A; headroom = A_max(M) − A.

Files (under perceptivity/data/):
  hpa_perceptivity_v03.npz   — R, A, C, headroom matrices (154 × 44)
  hpa_perceptivity_v03.json  — index/columns/A_max/cell_type_class metadata
"""

from functools import lru_cache
from pathlib import Path

import json
import numpy as np
import pandas as pd

_DATA_DIR = Path(__file__).parent / "data"
_NPZ = _DATA_DIR / "hpa_perceptivity_v03.npz"
_META = _DATA_DIR / "hpa_perceptivity_v03.json"


_NPZ_SCSCALED = _DATA_DIR / "hpa_perceptivity_v03_scscaled.npz"
_META_SCSCALED = _DATA_DIR / "hpa_perceptivity_v03_scscaled.json"

_OPINT = _DATA_DIR / "operation_intensity_v1.json"  # Factor-2 reference (v0.4)


@lru_cache(maxsize=4)
def load_hpa_perceptivity(mode: str = "pseudobulk"):
    """Load the bundled 154 × 44 HPA perceptivity reference.

    Parameters
    ----------
    mode : {"pseudobulk", "single_cell_scaled"}
        - "pseudobulk" (default, since v0.2.2): R, A computed directly as
          mean log1p(nCPM) over module genes per HPA cell type (pseudobulk units).
          Use when comparing to other pseudobulk data or for cell-type-level analyses.
        - "single_cell_scaled" (added in v0.2.3): HPA pseudobulk nCPM passed
          through scanpy normalize_total(1e4) + log1p before module scoring,
          matching the scale of per-cell scRNA-seq data after the same
          preprocessing. Use this mode when projecting single-cell tumor
          datasets per-cell into the eigenspace — guarantees that the per-cell
          cancer z-scores are commensurate with the HPA reference distribution.

    Returns
    -------
    dict
        R, A, C, headroom : DataFrame (154 × 44)
        A_max             : Series (44,)  per-module max A across HPA
        cell_type_class   : Series (154,) HPA "Cell type class" label
        mode              : str — which calibration was loaded
    """
    if mode not in ("pseudobulk", "single_cell_scaled"):
        raise ValueError(f"mode must be 'pseudobulk' or 'single_cell_scaled', got {mode!r}")

    if mode == "pseudobulk":
        npz_path, meta_path = _NPZ, _META
    else:
        npz_path, meta_path = _NPZ_SCSCALED, _META_SCSCALED

    if not npz_path.exists() or not meta_path.exists():
        if mode == "single_cell_scaled":
            raise FileNotFoundError(
                f"single-cell-scaled HPA reference not found at {npz_path}. "
                "Run scripts/12_build_hpa_scscaled.py to generate. "
                "Added in perceptome v0.2.3."
            )
        raise FileNotFoundError(
            f"HPA perceptivity reference not found at {npz_path}. "
            "Run scripts/02_build_hpa_perceptivity.py to generate."
        )

    arrs = np.load(npz_path, allow_pickle=False)
    with open(meta_path) as f:
        meta = json.load(f)

    cell_types = meta["cell_types"]
    modules = meta["modules"]
    A_max = pd.Series(meta["A_max"], index=modules, name="A_max")
    cell_type_class = pd.Series(meta["cell_type_class"], index=cell_types, name="cell_type_class")

    R = pd.DataFrame(arrs["R"], index=cell_types, columns=modules)
    A = pd.DataFrame(arrs["A"], index=cell_types, columns=modules)
    C = pd.DataFrame(arrs["C"], index=cell_types, columns=modules)
    headroom = pd.DataFrame(arrs["headroom"], index=cell_types, columns=modules)

    return {
        "R": R, "A": A, "C": C, "headroom": headroom,
        "A_max": A_max, "cell_type_class": cell_type_class,
        "mode": mode,
    }


@lru_cache(maxsize=1)
def load_operation_intensity():
    """Load the bundled Factor-2 reference (operation_intensity_v1).

    Calibrated from the substrate-series corpus (papers 4.3-4.8 + paper4 memory)
    by scripts/14_build_operation_intensity.py. Ordinal intensity bands per
    (operation, module) — NOT a validated point estimate (see `disclaimer`).

    Returns
    -------
    dict
        meta  : dict (schema_version, definition, bands, disclaimer, operations,
                cell_to_hpa_map)
        table : dict {(operation, module): record}, where record has
                intensity_band, direction, peak_effect, effect_metric, n_obs,
                calib_saturated, peak_condition, calib_cell.
    """
    if not _OPINT.exists():
        raise FileNotFoundError(
            f"operation_intensity reference not found at {_OPINT}. "
            "Run scripts/14_build_operation_intensity.py to generate. "
            "Added in perceptome v0.4 (Consolidation)."
        )
    with open(_OPINT) as f:
        payload = json.load(f)
    table = {(r["operation"], r["module"]): r for r in payload["table"]}
    meta = {k: payload[k] for k in
            ("schema_version", "definition", "bands", "disclaimer", "operations",
             "cell_to_hpa_map") if k in payload}
    return {"meta": meta, "table": table}


def hpa_capacity_floor(hi=4.5, lo=2.5):
    """Cell-type × module DataFrame of capacity-floor labels for the bundled HPA reference."""
    from .floor import capacity_floor
    ref = load_hpa_perceptivity()
    A = ref["A"]
    out = A.copy().astype(object)
    for c in A.columns:
        out[c] = [capacity_floor(v, hi=hi, lo=lo) for v in A[c].values]
    return out
