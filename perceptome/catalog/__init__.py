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
"""

from .modules import load_catalog, list_modules, get_genes, get_module_info, add_module
from .validate import validate_catalog

__all__ = [
    "load_catalog", "list_modules", "get_genes", "get_module_info",
    "add_module", "validate_catalog",
]
