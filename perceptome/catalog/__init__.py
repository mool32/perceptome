"""Module catalog — 44 regulatory modules (earlier called perceptual modules).

The canonical catalog is loaded once per process via load_catalog().
``core_genes`` are the 206 unique genes of the paper's Table S9 (219 gene slots); the four operational
criteria per module (Table S9b) are in ``module_criteria()``. Gene symbols are kept as published;
renamed ones are resolved by ``aliases_of`` / ``score_modules(resolve_aliases=True)`` (Note S-Q).
Each module entry exposes:
  - core_genes      sensor cascade + TF (4-5 genes; the scored panel / readiness measurement)
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
from .aliases import ALIASES, aliases_of, resolve_genes
from .criteria import module_criteria

__all__ = [
    "load_catalog", "list_modules", "get_genes", "get_module_info",
    "add_module", "validate_catalog",
    "ALIASES", "aliases_of", "resolve_genes", "module_criteria",
]
