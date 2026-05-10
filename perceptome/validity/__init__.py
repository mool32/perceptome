"""Validity controls for module-score perturbation analyses.

Three nulls + one positive control. Standard practice since Paper 4.5 v1.2 amendment
(random-200 ARTIFACT FAIL detection caught proliferation-baseline confound that
biased mean_raw scoring; switched to scanpy_corrected as the fix).

Contracts (from paper4.7 PRIMARY pre-reg):
  random_200       |log2FC| < 0.10        baseline shift artifact detector
  housekeeping     |log2FC| < 0.20        technical normalization sanity
  cell_cycle       |log2FC| > 0.30        positive control: perturbation reaches biology

A scorecard of PASS/FAIL/ARTIFACT per check is the standard output.
"""

from .nulls import (
    random_200_panel,
    housekeeping_panel,
    cell_cycle_panel,
    log2fc_perturbation_vs_control,
)
from .scorecard import validate_perturbation, scorecard

__all__ = [
    "random_200_panel", "housekeeping_panel", "cell_cycle_panel",
    "log2fc_perturbation_vs_control",
    "validate_perturbation", "scorecard",
]
