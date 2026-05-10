"""perceptome — cellular perception analysis toolkit (v0.2.0).

Three layers of analysis:

  1. GEOMETRY   Where does the cell live in module space?
                  pct.score_modules()    44 modules per cell
                  pct.project()          → 9-PC eigenspace (HPA-derived)
                  pct.compare_to_references()  cosines vs disease/aging/attractor/drugs

  2. PERCEPTIVITY   What can the cell DO? (capacity layer, NEW in v0.2)
                  pct.compute_perceptivity()   R, A, C, headroom, BS, GS, I, quadrant
                  pct.predict_engagement()     two-factor framework, capacity floor
                  pct.load_hpa_perceptivity()  154 × 44 reference matrices

  3. VALIDITY   Did the measurement reach biology, or hit an artifact?
                  pct.validate_perturbation()  random-200 / housekeeping / cell-cycle scorecard

Quick start:
    import scanpy as sc
    import perceptome as pct

    adata = sc.read_h5ad("your_data.h5ad")
    scores = pct.score_modules(adata)["scores"]              # 44 modules per cell
    coords = pct.project(scores)                             # 9-PC eigenspace
    perc = pct.compute_perceptivity(
        pct.score_readiness(adata), pct.score_activity(adata),
        cell_type=adata.obs["cell_type"], cell_class=adata.obs["tissue"],
    )
    refs = pct.compare_to_references(coords["coordinates"])  # disease/aging/attractor

See docs and examples/ for end-to-end pipelines.
"""

__version__ = "0.2.0"

from .catalog import load_catalog, list_modules, get_genes, get_module_info, add_module, validate_catalog
from .score import score_modules, score_readiness, score_activity
from .perceptivity import (
    compute_perceptivity, perceptivity_per_celltype, classify_quadrant,
    capacity_floor, predict_engagement,
    load_hpa_perceptivity, hpa_capacity_floor,
)
from .eigenspace import project, rebuild
from .compare import compare_conditions, divergence_score, compare_to_references, infrastructure_regime
from .reference import load_attractor_direction, attractor_alignment
from .drugs import drug_anchors, activity_layer_screen
from .validity import (
    random_200_panel, housekeeping_panel, cell_cycle_panel,
    log2fc_perturbation_vs_control,
    validate_perturbation, scorecard,
)
from .network import compute_network, compare_networks, module_heterogeneity
from .utils import cosine

__all__ = [
    "__version__",
    # catalog
    "load_catalog", "list_modules", "get_genes", "get_module_info", "add_module", "validate_catalog",
    # score
    "score_modules", "score_readiness", "score_activity",
    # perceptivity
    "compute_perceptivity", "perceptivity_per_celltype", "classify_quadrant",
    "capacity_floor", "predict_engagement",
    "load_hpa_perceptivity", "hpa_capacity_floor",
    # eigenspace
    "project", "rebuild",
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
    # utils
    "cosine",
]
