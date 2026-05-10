"""Reference perturbation vectors and direction loaders.

All bundled references live under reference/data/:
  - disease_vectors.json     AD/DKD/RA/IPF perturbation vectors (12-PC, ported from v0.1)
                             NOTE: still 12-PC v0.2 eigenspace; recompute on v0.3 (44-mod, 9-PC)
                             scheduled in cleanup phase.
  - aging_reference.json     inflammaging + collapse axes (12-PC, ported from v0.1)
  - attractor_v1.json        cancer capacity-direction reference (NEW v0.3, Paper 4.2 P3)
  - drug_vectors.csv         CMap/LINCS L1000 vectors (will be recomputed on 44-mod last)
"""

from .attractor import load_attractor_direction, attractor_alignment

__all__ = ["load_attractor_direction", "attractor_alignment"]
