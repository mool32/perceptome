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


@lru_cache(maxsize=2)
def load_hpa_perceptivity():
    """Load the bundled 154 × 44 HPA perceptivity reference.

    Returns
    -------
    dict
        R, A, C, headroom : DataFrame (154 × 44)
        A_max             : Series (44,)  per-module max A across HPA
        cell_type_class   : Series (154,) HPA "Cell type class" label
    """
    if not _NPZ.exists() or not _META.exists():
        raise FileNotFoundError(
            f"HPA perceptivity reference not found at {_NPZ}. "
            "Run scripts/02_build_hpa_perceptivity.py to generate."
        )

    arrs = np.load(_NPZ, allow_pickle=False)
    with open(_META) as f:
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
    }


def hpa_capacity_floor(hi=4.5, lo=2.5):
    """Cell-type × module DataFrame of capacity-floor labels for the bundled HPA reference."""
    from .floor import capacity_floor
    ref = load_hpa_perceptivity()
    A = ref["A"]
    out = A.copy().astype(object)
    for c in A.columns:
        out[c] = [capacity_floor(v, hi=hi, lo=lo) for v in A[c].values]
    return out
