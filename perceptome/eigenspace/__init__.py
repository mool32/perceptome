"""Reference eigenspace — the fixed 9-PC coordinate system built from 154 normal HPA cell types.

The 154 x 44 matrix of module scores (Human Protein Atlas, normal human cell types)
is z-scored per module and decomposed by correlation PCA; nine components are
retained by the Kaiser criterion (73.8% of variance; PC1-PC6 bootstrap-stable).
No cancer data enter any step, and the space is never refit: cells from any dataset
are projected into it (``project``) after standardization against a reference
calibration (``load_calibration``; pseudobulk bundled, single-cell-scaled for per-cell data).

For reproducibility of the first Paper 3 draft, the 43-module 12-PC eigenspace is
documented in docs/EIGENSPACE_EVOLUTION.md.
"""

from .project import project, project_shift
from .rebuild import rebuild
from .calibration import Calibration, load_calibration, hpa_reference_scores, read_hpa_single_cell_type

__all__ = [
    "project", "project_shift", "rebuild",
    "Calibration", "load_calibration", "hpa_reference_scores", "read_hpa_single_cell_type",
]
