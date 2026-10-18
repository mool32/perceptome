"""Convergence of malignant cells onto one region of normal cell-state space (Paper 3).

Spiro (2026), *Diverse cancers converge on one region of normal cell-state space in a
cancer-blind coordinate system* (Zenodo 10.5281/zenodo.20542130). The package implements the
paper's assay on any single-cell cohort:

  assign_cells / assign_to_reference   QC, score 44 modules, standardize, project, assign each cell
                                       to its nearest of 152 normal cell types (cosine)
  core_fraction / region_occupancy     occupancy and enrichment of the core / neighbourhood / region
  within_patient_occupancy             paired tumour-vs-matched-normal test (Wilcoxon + sign test)
  pooled_within_patient                the same pooled over datasets: mixed model, dataset random intercept
  occupant_contrast                    malignant vs normal occupants of the state, module by module (Fig. 6)
  contraction / contraction_null       mean pairwise cosine distance of cancer centroids vs their origins,
                                       and the matched-averaging null (Fig. 2, Note S-G)
  recurrence_ranking                   cross-cancer recurrence, depth and occupancy (Table S2)
  core_signature / target_distinctiveness / module_contributions
                                       what defines the core and how a cancer reproduces it
  load_convergence_reference           the paper's definitions, published results and SI tables

What the paper supports and what it does not (static convergence, not a proven attractor; the
region is robust, the cytotrophoblast anchor is catalog-conditional) is in docs/SCOPE.md.
Per-cell data need the single-cell-scaled calibration — see docs/CONVERGENCE.md.
"""

from .reference import (
    load_convergence_reference, published_within_patient, region_members, quarantined_cell_types,
)
from .qc import prepare_expression, detected_genes
from .assign import (
    assign_to_reference, assign_cells, in_region, core_fraction, enrichment, region_occupancy,
)
from .occupancy import within_patient_occupancy, recurrence_ranking
from .pooled import pooled_within_patient
from .contrast import occupant_contrast
from .cross_origin import (
    contraction, contraction_null, cosine_distances, mean_pairwise_cosine_distance,
)
from .signature import (
    core_profile, core_signature, target_distinctiveness, module_contributions,
    check_calibration_against_table_s7, S7_TOLERANCE, S7_KNOWN_DEVIATIONS,
)

__all__ = [
    "load_convergence_reference", "published_within_patient", "region_members", "quarantined_cell_types",
    "prepare_expression", "detected_genes",
    "assign_to_reference", "assign_cells", "in_region", "core_fraction", "enrichment", "region_occupancy",
    "within_patient_occupancy", "recurrence_ranking", "pooled_within_patient", "occupant_contrast",
    "contraction", "contraction_null", "cosine_distances", "mean_pairwise_cosine_distance",
    "core_profile", "core_signature", "target_distinctiveness", "module_contributions",
    "check_calibration_against_table_s7", "S7_TOLERANCE", "S7_KNOWN_DEVIATIONS",
]
