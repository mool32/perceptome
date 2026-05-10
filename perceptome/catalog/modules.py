"""Catalog access — load, list, query gene sets per module."""

import json
from functools import lru_cache
from pathlib import Path

_CATALOG_FILE = Path(__file__).parent / "data" / "modules_v03.json"


@lru_cache(maxsize=2)
def load_catalog(path=None):
    """Load module catalog.

    Parameters
    ----------
    path : str | Path | None
        Override path. Default = bundled v0.3 catalog (44 modules).

    Returns
    -------
    dict
        Top-level keys: 'version', 'n_modules', 'modules', 'release_notes'.
        modules is a dict {module_name: {core_genes, activity_genes, ...}}.
    """
    p = Path(path) if path else _CATALOG_FILE
    with open(p) as f:
        return json.load(f)


def list_modules(catalog=None, pan_cellular_only=False):
    """Return sorted module names.

    Parameters
    ----------
    catalog : dict | None
    pan_cellular_only : bool
        If True, exclude modules with pan_cellular=False (none in v0.3,
        but reserved for v0.4 tissue-specific extensions).
    """
    if catalog is None:
        catalog = load_catalog()
    names = catalog["modules"].keys()
    if pan_cellular_only:
        names = [n for n in names if catalog["modules"][n].get("pan_cellular", True)]
    return sorted(names)


def get_genes(module, gene_set="core", catalog=None):
    """Return gene list for a module.

    Parameters
    ----------
    module : str
    gene_set : 'core' | 'activity' | 'sensor' | 'cascade' | 'tf' | 'feedback' | 'full'
        - core      readiness genes (machinery)
        - activity  TF target genes (engagement)  — recommended for A in perceptivity
        - full      union of all gene categories
    catalog : dict | None
    """
    if catalog is None:
        catalog = load_catalog()
    mod = catalog["modules"].get(module)
    if mod is None:
        raise KeyError(f"Module '{module}' not found. Available: {list_modules(catalog)}")

    if gene_set == "core":
        return list(mod["core_genes"])
    if gene_set == "activity":
        return list(mod.get("activity_genes", mod["core_genes"]))
    if gene_set == "sensor":
        return list(mod.get("sensor_genes", []))
    if gene_set == "cascade":
        return list(mod.get("cascade_genes", []))
    if gene_set == "tf":
        return list(mod.get("tf_genes", mod.get("primary_tf", [])))
    if gene_set == "feedback":
        return list(mod.get("feedback_genes", []))
    if gene_set == "full":
        bag = set(mod["core_genes"])
        for k in ("sensor_genes", "cascade_genes", "tf_genes", "feedback_genes", "activity_genes"):
            bag.update(mod.get(k, []))
        return sorted(bag)

    raise ValueError(
        f"gene_set must be one of core|activity|sensor|cascade|tf|feedback|full, got {gene_set!r}"
    )


def get_module_info(module, catalog=None):
    """Return full info dict for a module (copy)."""
    if catalog is None:
        catalog = load_catalog()
    mod = catalog["modules"].get(module)
    if mod is None:
        raise KeyError(f"Module '{module}' not found.")
    return dict(mod)


def add_module(name, definition, catalog_path=None):
    """Append a module to a catalog file and re-save.

    Use for ad-hoc extensions in user code; for the canonical bundle, edit
    the catalog build script and re-run instead.
    """
    p = Path(catalog_path) if catalog_path else _CATALOG_FILE
    load_catalog.cache_clear()
    with open(p) as f:
        catalog = json.load(f)
    catalog["modules"][name] = definition
    catalog["n_modules"] = len(catalog["modules"])
    with open(p, "w") as f:
        json.dump(catalog, f, indent=2, ensure_ascii=False)
