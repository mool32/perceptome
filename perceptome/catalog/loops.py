"""Loop layer (v0.4) — negative-feedback structure of each module.

Sits next to the v0.3 module catalog without changing it. Built by
scripts/12_build_loops_v04.py; see docs/LOOPS.md for field definitions.

Status: curated, not yet validated.
"""

import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

from .modules import load_catalog

_LOOPS_FILE = Path(__file__).parent / "data" / "loops_v04.json"
_TIERS = ("A", "B", "C", "D", "E")
_PERTURBATION_CLASSES = ("break", "input_down", "input_up")


@lru_cache(maxsize=2)
def load_loops(path=None):
    """Load the loop layer. Default = bundled loops_v04.json (44 modules)."""
    p = Path(path) if path else _LOOPS_FILE
    with open(p) as f:
        return json.load(f)


def list_loops(tiers=("A", "B", "C", "D"), loops=None):
    """Return sorted module names whose loop tier is in `tiers`.

    Default excludes tier E (no clean transcriptional loop).
    """
    if loops is None:
        loops = load_loops()
    if isinstance(tiers, str):
        tiers = (tiers,)
    return sorted(n for n, l in loops["loops"].items() if l["tier"] in tiers)


def get_loop(module, loops=None):
    """Return the loop entry for a module (copy)."""
    if loops is None:
        loops = load_loops()
    entry = loops["loops"].get(module)
    if entry is None:
        raise KeyError(f"Module '{module}' not in loop layer. Available: {sorted(loops['loops'])}")
    return json.loads(json.dumps(entry))


def loop_gene_sets(module, loops=None):
    """Return {'feedback': [...], 'targets': [...]} — disjoint gene sets for scoring."""
    entry = get_loop(module, loops)
    return {"feedback": list(entry["feedback_core"]), "targets": list(entry["loop_targets"])}


def loop_output_genes(module, loops=None):
    """Genes that read out the module's output.

    loop_targets, plus feedback_core when the feedback IS the output
    (chaperone titration, clock, ligand catabolism): there the feedback genes
    are the module's main induced products, so excluding them would discard
    the readout.
    """
    entry = get_loop(module, loops)
    out = list(entry["loop_targets"])
    if entry["feedback_is_output"]:
        out = list(entry["feedback_core"]) + [g for g in out if g not in entry["feedback_core"]]
    return out


def loop_overlaps(loops=None, field="loop_targets"):
    """Genes that appear in `field` of more than one module: {gene: [modules]}.

    Shared genes are biology, not errors, but a metric that claims module
    specificity must account for them.
    """
    if loops is None:
        loops = load_loops()
    where = defaultdict(list)
    for name, entry in loops["loops"].items():
        for g in entry[field]:
            where[g].append(name)
    return {g: sorted(m) for g, m in sorted(where.items()) if len(m) > 1}


def validate_loops(loops=None, catalog=None, strict=False):
    """Return list of (severity, module, message) tuples. Empty list = clean.

    Errors: unknown module, bad tier, feedback/target overlap, alias symbols,
    missing perturbation classes. Warnings: tier A/B loops with < 3 targets or
    without any readable feedback and no break perturbation.
    """
    if loops is None:
        loops = load_loops()
    if catalog is None:
        catalog = load_catalog()
    aliases = set(loops.get("aliases", {}))
    issues = []

    missing = set(catalog["modules"]) - set(loops["loops"])
    for name in sorted(missing):
        issues.append(("error", name, "module in catalog but not in loop layer"))

    for name, e in loops["loops"].items():
        if name not in catalog["modules"]:
            issues.append(("error", name, "module not in base catalog"))
        if e["tier"] not in _TIERS:
            issues.append(("error", name, f"tier {e['tier']!r} not in {_TIERS}"))
        fb, tg = set(e["feedback_core"]), set(e["loop_targets"])
        if fb & tg:
            issues.append(("error", name, f"feedback_core overlaps loop_targets: {sorted(fb & tg)}"))
        pert = e.get("perturbations", {})
        if set(pert) != set(_PERTURBATION_CLASSES):
            issues.append(("error", name, f"perturbations must have keys {_PERTURBATION_CLASSES}"))
        genes = fb | tg | {g for v in pert.values() for g in v}
        if genes & aliases:
            issues.append(("error", name, f"alias symbols used: {sorted(genes & aliases)}"))
        if e["tier"] in ("A", "B"):
            if len(tg) < 3:
                issues.append(("warn", name, f"tier {e['tier']} loop with {len(tg)} targets (< 3)"))
            if not fb and not pert.get("break") and not pert.get("input_up"):
                issues.append(("warn", name, "no readable feedback and no break/input_up perturbation"))

    if strict and any(s == "error" for s, _, _ in issues):
        raise ValueError(f"Loop layer validation failed: {issues}")
    return issues
