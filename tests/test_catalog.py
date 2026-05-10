"""Catalog tests — 44 modules, NPAS4 present, validation clean."""

import perceptome as pct


def test_version_string():
    # Version follows semver; check format rather than exact value to avoid
    # churning this test on every patch release.
    parts = pct.__version__.split(".")
    assert len(parts) >= 2 and all(p.isdigit() or "-" in p for p in parts), pct.__version__
    assert pct.__version__.startswith("0.2.")


def test_catalog_has_44_modules():
    cat = pct.load_catalog()
    assert cat["n_modules"] == 44
    assert len(cat["modules"]) == 44
    assert cat["version"] == "0.3"


def test_npas4_in_catalog():
    cat = pct.load_catalog()
    assert "NPAS4" in cat["modules"]
    npas4 = cat["modules"]["NPAS4"]
    assert "NPAS4" in npas4["core_genes"]
    assert "BDNF" in npas4["activity_genes"]
    assert npas4["dissociation_risk"] == "HIGH"
    assert npas4["category"] == "A_exteroceptive"
    assert "neuron" in npas4["tissue_bias"]


def test_get_genes_all_sets():
    for s in ("core", "activity", "sensor", "cascade", "tf", "feedback", "full"):
        out = pct.get_genes("NF-κB", gene_set=s)
        assert isinstance(out, list)
    full = set(pct.get_genes("NPAS4", "full"))
    core = set(pct.get_genes("NPAS4", "core"))
    activity = set(pct.get_genes("NPAS4", "activity"))
    assert core.issubset(full) and activity.issubset(full)


def test_validate_catalog_clean():
    issues = pct.validate_catalog()
    errors = [i for i in issues if i[0] == "error"]
    warnings = [i for i in issues if i[0] == "warn"]
    assert errors == []
    assert warnings == []


def test_get_genes_unknown_module():
    import pytest
    with pytest.raises(KeyError):
        pct.get_genes("NonExistentModule", "core")


def test_list_modules_pan_cellular_filter():
    all_mods = pct.list_modules()
    pan = pct.list_modules(pan_cellular_only=True)
    assert len(all_mods) == 44
    assert len(pan) == 44  # all v0.3 modules are pan-cellular currently
