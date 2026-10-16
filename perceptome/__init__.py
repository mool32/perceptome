"""perceptome — the cancer-blind coordinate system of cell state (44 regulatory modules).

Implements Spiro (2026), "Diverse cancers converge on one region of normal cell-state space in a
cancer-blind coordinate system" (Zenodo 10.5281/zenodo.20542130), and the capacity / validity /
drug layers of the companion papers. Four layers of analysis:

  1. GEOMETRY      Where does the cell live among normal cell states?
                     pct.score_modules()      44 regulatory modules per cell (mean_raw)
                     pct.project()            standardize against a calibration, project into the
                                              fixed 9-PC space built from 154 normal HPA cell types
                     pct.compare_to_references()  cosines vs legacy cancer / disease / aging directions

  2. CONVERGENCE   Do malignant cells sit in one region of that space, and more so than the
                   same patient's normal cells?  (pct.convergence, Paper 3)
                     assign_cells()           QC -> score -> project -> nearest of 152 normal cell types
                     within_patient_occupancy()   paired tumour vs matched-normal test
                     load_convergence_reference() the paper's definitions, results and SI tables

  3. PERCEPTIVITY  What can the cell DO? (capacity layer)
                     pct.compute_perceptivity()   R, A, C, headroom, BS, GS, I, quadrant
                     pct.predict_engagement()     two-factor framework, capacity floor
                     pct.load_hpa_perceptivity()  154 x 44 reference matrices

  4. VALIDITY      Did the measurement reach biology, or hit an artifact?
                     pct.validate_perturbation()  random-200 / housekeeping / cell-cycle scorecard

Quick start (per-cell data need the single_cell_scaled calibration; see docs/CONVERGENCE.md):
    import scanpy as sc
    import perceptome as pct

    adata = sc.read_h5ad("cohort.h5ad")                      # tumour + matched normal, HGNC symbols
    res = pct.convergence.assign_cells(adata, calibration="single_cell_scaled")
    test = pct.within_patient_occupancy(
        res["assignment"], adata.obs["patient"], adata.obs["arm"], malignant="tumour", normal="normal")

    scores = pct.score_modules(adata)["scores"]              # 44 modules per cell
    coords = pct.project(scores)["coordinates"]              # standardized 9-PC coordinates

See docs/ and examples/ for end-to-end pipelines, docs/SCOPE.md for what the framework does and
does not support, and docs/PAPER_MAP.md for how the paper maps onto the code.
"""

__version__ = "0.2.2"

from .catalog import (
    load_catalog, list_modules, get_genes, get_module_info, add_module, validate_catalog,
    aliases_of, resolve_genes, module_criteria,
)
from .score import score_modules, score_readiness, score_activity
from .perceptivity import (
    compute_perceptivity, perceptivity_per_celltype, classify_quadrant,
    capacity_floor, predict_engagement,
    load_hpa_perceptivity, hpa_capacity_floor,
)
from .eigenspace import (
    project, project_shift, rebuild, Calibration, load_calibration, hpa_reference_scores,
    read_hpa_single_cell_type,
)
from .compare import compare_conditions, divergence_score, compare_to_references, infrastructure_regime
from .reference import load_attractor_direction, attractor_alignment
from .drugs import drug_anchors, activity_layer_screen
from .validity import (
    random_200_panel, housekeeping_panel, cell_cycle_panel,
    log2fc_perturbation_vs_control,
    validate_perturbation, scorecard,
)
from .network import compute_network, compare_networks, module_heterogeneity
from . import convergence
from .convergence import (
    assign_cells, assign_to_reference, core_fraction, region_occupancy,
    within_patient_occupancy, recurrence_ranking, load_convergence_reference,
)
from .utils import cosine

__all__ = [
    "__version__",
    # catalog
    "load_catalog", "list_modules", "get_genes", "get_module_info", "add_module", "validate_catalog",
    "aliases_of", "resolve_genes", "module_criteria",
    # score
    "score_modules", "score_readiness", "score_activity",
    # perceptivity
    "compute_perceptivity", "perceptivity_per_celltype", "classify_quadrant",
    "capacity_floor", "predict_engagement",
    "load_hpa_perceptivity", "hpa_capacity_floor",
    # eigenspace
    "project", "project_shift", "rebuild",
    "Calibration", "load_calibration", "hpa_reference_scores", "read_hpa_single_cell_type",
    # compare
    "compare_conditions", "divergence_score", "compare_to_references", "infrastructure_regime",
    # reference
    "load_attractor_direction", "attractor_alignment",
    # drugs (narrow validated scope — see Paper 4.1)
    "drug_anchors", "activity_layer_screen",
    # validity
    "random_200_panel", "housekeeping_panel", "cell_cycle_panel",
    "log2fc_perturbation_vs_control",
    "validate_perturbation", "scorecard",
    # network
    "compute_network", "compare_networks", "module_heterogeneity",
    # convergence (Paper 3)
    "convergence", "assign_cells", "assign_to_reference", "core_fraction", "region_occupancy",
    "within_patient_occupancy", "recurrence_ranking", "load_convergence_reference",
    # utils
    "cosine",
]
