"""Perceptivity layer — capacity / engagement / specialization metrics.

Companion to v1.4 metric design (perceptivity_metric/METRIC_DESIGN_v1_4_amendment.md).

Per (cell type, module):
  R          readiness        mean log expression of core_genes
  A          activity         mean log expression of activity_genes
  C          capacity         R − A  (signed; positive = ramp room)
  headroom   distance to max  A_max(M) − A
  I          intensity        per-engaged-module mean A

Per cell type (aggregated across engaged modules with A > A_thr):
  BS         within-class specialization breadth (soft)
  GS         global specialization breadth (soft)
  spec_quadrant  Q1 unique | Q2 locally homogeneous | Q3 globally common | Q4 generalist

Capacity-floor predictor (paper4.5 + 4.6 + 4.7 + 4.8 closed):
  A_baseline > 4.5  ⇒  upward saturation (cannot ramp UP further). Downward suppression
                       by specific signaling (e.g. atRA → UPR) IS possible — predictor
                       is upward-asymmetric.
  A_baseline < 2.5  ⇒  capacious (ramp possible; magnitude operation-determined)
  2.5 ≤ A_baseline ≤ 4.5  ⇒  intermediate

Two-factor framework:
  engagement(cell, op, module) = capacity(cell, module) × operation_intensity(op, module)
  Tool covers Factor 1 only (HPA-derivable).
"""

from .compute import compute_perceptivity, perceptivity_per_celltype, classify_quadrant
from .floor import capacity_floor, predict_engagement
from .reference import load_hpa_perceptivity, hpa_capacity_floor

__all__ = [
    "compute_perceptivity", "perceptivity_per_celltype", "classify_quadrant",
    "capacity_floor", "predict_engagement",
    "load_hpa_perceptivity", "hpa_capacity_floor",
]
