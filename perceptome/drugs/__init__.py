"""Drug perturbation analysis — narrow validated scope.

WARNING — read before using.

Paper 4.1 (perceptome on drugs) closed 2026-05-09 with 6 surviving validated
findings and 11 pre-registered falsifications (see PAPER4_1_CANONICAL_RESULTS.md
in the source tree, and killed_hypotheses.csv).

Operations FALSIFIED — do NOT do these in your work:
  - drug-disease cosine matching in eigenspace (3 formulations falsified)
  - drug-class mechanism deconvolution via panel geometry
  - readiness-layer timescale rescue (24h vs 6h)
  - TF-autoregulation as 1st-class layer
  - equilibrium-instrument linear-response framing
  - "clean perceptome signature" predicts FDA approval beyond selectivity
  - snapshot perceptome blind to pulsatile dynamics

Operation VALIDATED — narrow scope, this module supports it:
  Activity-layer scoring of TF target panels for 9 specific (class, module)
  anchors. Procedure: per-class observed activity z vs background null of
  non-panel drugs, BH-FDR q<0.10 (Block 5 v1.2, locked 2026-05-07).
  Discovery + holdout PASS for 2 cleanest anchors (MEKi/ERK, Proteasomei/UPR-PERK).

Tool API:
  pct.drug_anchors()                  9-row reference table with validation chain
  pct.activity_layer_screen(adata,    Per-class background-null test on user data
                            ...)
"""

from .anchors import drug_anchors, ANCHORS
from .screen import activity_layer_screen

__all__ = ["drug_anchors", "ANCHORS", "activity_layer_screen"]
