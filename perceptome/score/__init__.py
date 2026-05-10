"""Module activity scoring from expression data.

Two complementary measurements per (cell, module):
  - readiness  R = mean expression of core_genes (machinery present)
  - activity   A = mean expression of activity_genes (machinery currently engaged)

Higher-level perceptivity quantities (R, A, C, headroom, BS, GS, I) are computed
in pct.perceptivity from these primitives.
"""

from .methods import score_modules, score_readiness, score_activity

__all__ = ["score_modules", "score_readiness", "score_activity"]
