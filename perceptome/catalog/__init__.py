"""Module catalog — 44 perceptual signaling modules.

The canonical catalog is loaded once per process via load_catalog().
Each module entry exposes:
  - core_genes      sensor cascade + TF (≈3-5 genes; readiness measurement)
  - activity_genes  TF target genes (≈5-15 genes; engagement measurement)
  - primary_tf      the dominant TF(s) for the module
  - sensor_genes / cascade_genes / tf_genes / feedback_genes  optional structural detail
  - category        one of A_exteroceptive | A_interoceptive | B_interoceptive
                          | nuclear_receptor | infrastructure
  - dissociation_risk      LOW | MEDIUM | HIGH (R vs A divergence likelihood)
  - dissociation_note      explanation when MEDIUM/HIGH
  - pan_cellular           True if module operates across all cell types
  - tissue_bias            list of cell-class hints (does NOT exclude scoring)
  - mii                    Module Importance Index (or None if not yet computed)

Loop layer (v0.4, separate file, v0.3 fields untouched): per-module negative-
feedback structure — feedback_core, loop_targets (disjoint), perturbation
classes and a testability tier. See load_loops() and docs/LOOPS.md.
"""

from .modules import load_catalog, list_modules, get_genes, get_module_info, add_module
from .validate import validate_catalog
from .loops import (
    load_loops, list_loops, get_loop, loop_gene_sets, loop_overlaps, validate_loops,
)

__all__ = [
    "load_catalog", "list_modules", "get_genes", "get_module_info",
    "add_module", "validate_catalog",
    "load_loops", "list_loops", "get_loop", "loop_gene_sets", "loop_overlaps", "validate_loops",
]
