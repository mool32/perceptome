"""What defines the convergence core in module space, and how cancers reproduce it.

The core is defined from normal cells alone. In the standardized module space (reference mean 0,
SD 1) the two-component signature is the unit-normalized mean of the cytotrophoblast and
megakaryocyte profiles; a cancer's per-module contribution to its proximity to the core is the
elementwise term of the cosine between its unit-normalized standardized centroid and that
signature, so the 44 contributions sum to the cosine (Table S4; mean 0.711 over the 25 cancers).
"""

import numpy as np
import pandas as pd

from ..eigenspace.calibration import Calibration, load_calibration
from .reference import load_convergence_reference, region_members


def _calibration(calibration):
    return calibration if isinstance(calibration, Calibration) else load_calibration(calibration)


def core_profile(calibration="pseudobulk", region="core"):
    """Standardized module profile (z per module) of each member of ``region`` and their mean.

    Returns a DataFrame (modules x [members..., 'mean']). Table S7 prints this for the
    megakaryocyte and cytotrophoblast under the paper's single-cell-scaled calibration.
    """
    cal = _calibration(calibration)
    members = region_members(region)
    z = cal.standardize(cal.reference_scores.loc[members]).T
    z["mean"] = z.mean(axis=1)
    return z


def core_signature(calibration="pseudobulk", region="core"):
    """Unit-normalized mean standardized profile of the region (44-vector, Series)."""
    mean = core_profile(calibration, region)["mean"]
    return mean / np.linalg.norm(mean.values)


def module_contributions(scores, calibration="pseudobulk", region="core"):
    """Per-module contribution of each row of module scores to its cosine with the core signature.

    ``scores`` is one or several module-score vectors (e.g. a cancer's malignant centroid).
    Rows are standardized against the calibration, unit-normalized, and multiplied elementwise by
    the signature; each row of the result sums to the cosine.
    """
    cal = _calibration(calibration)
    single = isinstance(scores, pd.Series)
    X = scores.to_frame().T if single else scores
    z = cal.standardize(X)
    u = z.div(np.linalg.norm(z.values, axis=1), axis=0)
    out = u.mul(core_signature(cal, region), axis=1)
    return out.iloc[0] if single else out


def check_calibration_against_table_s7(calibration):
    """Compare a calibration's core z-profiles with the paper's Table S7.

    Table S7 prints, for each of the 44 modules, the z-score of the megakaryocyte and of the
    cytotrophoblast against the 154 normal cell types under the paper's single-cell-scaled
    calibration (two decimals). A calibration built to reproduce the paper should match them to
    rounding; use this after building ``single_cell_scaled`` from the HPA table.

    Returns
    -------
    dict
        megakaryocyte_r, cytotrophoblast_r   Pearson r over the 44 modules
        max_abs_diff, mean_abs_diff          over both cells and all modules
        matches                              True if max_abs_diff <= 0.011 (rounding of the print)
    """
    cal = _calibration(calibration)
    s7 = load_convergence_reference()["tables"]["S7"]
    z = cal.standardize(cal.reference_scores.loc[["megakaryocytes", "cytotrophoblasts"]]).T.reindex(s7.index)
    mk, ct = z["megakaryocytes"].values, z["cytotrophoblasts"].values
    d = np.abs(np.concatenate([mk - s7["megakaryocyte_z"].values, ct - s7["cytotrophoblast_z"].values]))
    return {
        "megakaryocyte_r": float(np.corrcoef(mk, s7["megakaryocyte_z"].values)[0, 1]),
        "cytotrophoblast_r": float(np.corrcoef(ct, s7["cytotrophoblast_z"].values)[0, 1]),
        "max_abs_diff": float(d.max()), "mean_abs_diff": float(d.mean()),
        "matches": bool(d.max() <= 0.011),
    }
