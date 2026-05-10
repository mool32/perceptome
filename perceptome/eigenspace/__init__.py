"""Perceptome eigenspace — 12-PC reference space derived from 154 × 44 HPA matrix.

For paper3 reproducibility, the original 12-PC v0.2 (43-module) eigenspace is
preserved in tool/perceptome/perceptome/data/. The v0.3 eigenspace adds NPAS4
as the 44th module and reruns the eigendecomposition on 154 × 44.
"""

from .project import project
from .rebuild import rebuild

__all__ = ["project", "rebuild"]
