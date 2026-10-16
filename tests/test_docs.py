"""Documentation integrity: relative links resolve; notebooks are valid, error-free, and the offline recipe ships with outputs."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = [ROOT / n for n in ("README.md", "CHANGELOG.md", "ROADMAP.md", "CONTRIBUTING.md")] + sorted((ROOT / "docs").glob("*.md"))
NOTEBOOKS = sorted((ROOT / "examples").glob("*.ipynb"))
LINK = re.compile(r"\]\(([^)\s]+)\)")


@pytest.mark.parametrize("doc", DOCS, ids=lambda p: p.name)
def test_relative_links_resolve(doc):
    targets = [m.group(1) for m in LINK.finditer(doc.read_text(encoding="utf-8"))]
    local = [t for t in targets if not t.startswith(("http://", "https://", "mailto:", "#"))]
    missing = [t for t in local if t.split("#")[0] and not (doc.parent / t.split("#")[0]).exists()]
    assert not missing, f"{doc.name} links to files that do not exist: {missing}"


@pytest.mark.parametrize("notebook", NOTEBOOKS, ids=lambda p: p.name)
def test_notebooks_are_valid_and_carry_no_error_output(notebook):
    nbformat = pytest.importorskip("nbformat")
    nb = nbformat.read(notebook, as_version=4)
    nbformat.validate(nb)
    errors = [(i, o.get("ename")) for i, c in enumerate(nb.cells) if c.cell_type == "code"
              for o in c.get("outputs", []) if o.get("output_type") == "error"]
    assert not errors, f"{notebook.name} stores error outputs: {errors}"


def test_the_offline_convergence_recipe_ships_with_outputs():
    """It needs no download, so every code cell is stored executed (the PBMC3K notebooks cannot promise that)."""
    nbformat = pytest.importorskip("nbformat")
    nb = nbformat.read(ROOT / "examples" / "recipe_convergence_occupancy.ipynb", as_version=4)
    code = [c for c in nb.cells if c.cell_type == "code"]
    assert code and all(c.get("execution_count") and c.get("outputs") for c in code)


def test_every_documented_builder_script_exists():
    text = "\n".join(p.read_text(encoding="utf-8") for p in DOCS)
    scripts = set(re.findall(r"scripts/(\d\d_[A-Za-z0-9_]+\.py)", text))
    missing = sorted(s for s in scripts if not (ROOT / "scripts" / s).exists())
    assert scripts and not missing, f"documents refer to scripts that do not exist: {missing}"
