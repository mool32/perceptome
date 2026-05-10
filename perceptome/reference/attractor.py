"""Cancer capacity-direction attractor reference (Paper 4.2 P3 PASS).

The attractor direction is a vector in (a) the 44-module per-module Δ space
and (b) the v0.3 12-PC eigenspace. Use it to test whether a (tumor − normal)
shift vector aligns with the observed cancer-transformation direction.

P3 PASS criterion (locked): cosine(shift, attractor_direction) > +0.20
in ≥ 4 of 6 cell types per cohort. Sun 2021 HCC paired cohort PASS.
"""

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ..utils import cosine

_REF_FILE = Path(__file__).parent / "data" / "attractor_v1.json"


@lru_cache(maxsize=1)
def load_attractor_direction():
    """Load the bundled attractor reference.

    Returns
    -------
    dict
        version
        attractor_cluster_cells          list[str]  8 HPA cell types
        attractor_direction_modules      Series (modules → Δ A_baseline)
        attractor_direction_eigenspace   ndarray (n_pcs,)
        eigenspace_PC_names              list[str]
        P3_per_celltype                  dict (Sun 2021 detail)
        source                           pre-reg sha + closing memo paths
    """
    with open(_REF_FILE) as f:
        d = json.load(f)
    d["attractor_direction_modules"] = pd.Series(
        d["attractor_direction_modules"], name="attractor_delta"
    )
    d["attractor_direction_eigenspace"] = np.array(d["attractor_direction_eigenspace"])
    return d


def attractor_alignment(shift_vector, mode="eigenspace"):
    """Cosine of a shift vector with the attractor direction.

    Parameters
    ----------
    shift_vector : Series | ndarray
        Tumor − normal Δ. If mode='modules', indexed by module name (44 modules).
        If mode='eigenspace', a vector in the 9-PC space.
    mode : 'eigenspace' | 'modules'

    Returns
    -------
    dict
        cosine                  float in [-1, +1]
        passes_p3_threshold     bool   cosine > +0.20
        threshold               0.20  (locked from Paper 4.2 pre-reg)
    """
    ref = load_attractor_direction()
    if mode == "eigenspace":
        attr = ref["attractor_direction_eigenspace"]
        if hasattr(shift_vector, "values"):
            shift = np.asarray(shift_vector.values, dtype=float)
        else:
            shift = np.asarray(shift_vector, dtype=float)
        attr = attr[: len(shift)]
        shift = shift[: len(attr)]
    elif mode == "modules":
        attr_s = ref["attractor_direction_modules"]
        if hasattr(shift_vector, "loc"):
            common = attr_s.index.intersection(shift_vector.index)
            attr = attr_s.loc[common].values.astype(float)
            shift = shift_vector.loc[common].values.astype(float)
        else:
            attr = attr_s.values
            shift = np.asarray(shift_vector, dtype=float)
    else:
        raise ValueError(f"mode must be 'eigenspace' or 'modules', got {mode!r}")

    cos = cosine(shift, attr)
    return {"cosine": float(cos), "passes_p3_threshold": cos > 0.20, "threshold": 0.20}
