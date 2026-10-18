"""Reference perturbation vectors and direction loaders.

All bundled references live under reference/data/ and carry their eigenspace vectors in the
standardized reference coordinates (``pct.project`` / ``pct.project_shift``; see
scripts/12_restandardize_reference_vectors.py):
  - disease_vectors.json     AD/DKD/RA/IPF perturbation vectors (9-PC, recomputed from raw scRNA-seq)
  - aging_reference.json     inflammaging + collapse axes (9-PC)
  - attractor_v1.json        legacy cancer capacity-direction reference (Paper 4.2 P3). Paper 3's
                             convergence result (25 cancers, nearest-cell assignment) lives in
                             perceptome.convergence and supersedes the 8-cell cluster.
"""

from .attractor import load_attractor_direction, attractor_alignment

__all__ = ["load_attractor_direction", "attractor_alignment"]
