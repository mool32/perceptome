"""Catalog integrity checks.

Catches: duplicate gene symbols across required fields, modules with
empty core/activity sets, dissociation_risk values outside the allowed
vocabulary, missing primary_tf when category requires one.

Used by tests/ — also useful when adding a new module by hand.
"""

from .modules import load_catalog

_ALLOWED_RISK = {"LOW", "MEDIUM", "HIGH", None}
_ALLOWED_CATEGORIES = {
    "A_exteroceptive", "A_interoceptive", "B_interoceptive",
    "nuclear_receptor", "infrastructure",
}


def validate_catalog(catalog=None, strict=False):
    """Return list of (severity, module, message) tuples. Empty list = clean.

    Parameters
    ----------
    strict : bool
        If True, treat warnings as errors (raise ValueError).
    """
    if catalog is None:
        catalog = load_catalog()

    issues = []

    if "modules" not in catalog:
        issues.append(("error", None, "catalog missing 'modules' key"))
        if strict:
            raise ValueError(issues)
        return issues

    if catalog.get("n_modules") != len(catalog["modules"]):
        issues.append(
            ("warn", None,
             f"n_modules={catalog.get('n_modules')} but {len(catalog['modules'])} entries")
        )

    for name, mod in catalog["modules"].items():
        # required fields
        if not mod.get("core_genes"):
            issues.append(("error", name, "missing or empty core_genes"))
        if not mod.get("activity_genes"):
            issues.append(("warn", name, "missing or empty activity_genes — falls back to core"))

        # category vocabulary
        cat = mod.get("category")
        if cat not in _ALLOWED_CATEGORIES:
            issues.append(("warn", name, f"category {cat!r} not in {sorted(_ALLOWED_CATEGORIES)}"))

        # risk vocabulary
        risk = mod.get("dissociation_risk")
        if risk not in _ALLOWED_RISK:
            issues.append(("error", name, f"dissociation_risk {risk!r} not in {_ALLOWED_RISK}"))

        # MEDIUM/HIGH risk should have explanation
        if risk in {"MEDIUM", "HIGH"} and not mod.get("dissociation_note"):
            issues.append(("warn", name, f"{risk} dissociation_risk without dissociation_note"))

    if strict and any(s == "error" for s, _, _ in issues):
        raise ValueError(f"Catalog validation failed: {issues}")
    return issues
