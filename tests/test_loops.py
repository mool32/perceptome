"""Loop layer (v0.4) — structure and integrity; v0.3 catalog untouched."""

import pytest

import perceptome as pct


def test_loop_layer_covers_catalog():
    loops = pct.load_loops()
    assert loops["base_catalog"] == "0.3"
    assert set(loops["loops"]) == set(pct.list_modules())


def test_validate_loops_clean():
    assert pct.validate_loops() == []


def test_feedback_and_targets_disjoint_everywhere():
    for name in pct.list_modules():
        sets = pct.loop_gene_sets(name)
        assert not set(sets["feedback"]) & set(sets["targets"]), name


def test_feedback_cores_not_shared_between_modules():
    assert pct.loop_overlaps(field="feedback_core") == {}


def test_tier_a_loops():
    tier_a = pct.list_loops("A")
    assert tier_a == ["ERK/MAPK", "HIF", "HSF1", "JAK-STAT", "NRF2", "SREBP", "UPR-ATF6", "UPR-PERK"]
    for name in tier_a:
        e = pct.get_loop(name)
        assert e["perturbations"]["break"] or e["perturbations"]["input_up"], name
        assert len(e["loop_targets"]) >= 3, name


def test_unreadable_feedback_has_no_feedback_genes():
    for name in pct.list_modules():
        e = pct.get_loop(name)
        if not e["feedback_readable"] and e["tier"] != "E":
            assert e["feedback_core"] == [], name


def test_no_alias_symbols_in_loop_layer():
    loops = pct.load_loops()
    aliases = set(loops["aliases"])
    for name, e in loops["loops"].items():
        genes = set(e["feedback_core"]) | set(e["loop_targets"])
        genes |= {g for v in e["perturbations"].values() for g in v}
        assert not genes & aliases, name


def test_validate_loops_catches_overlap():
    loops = pct.load_loops()
    broken = {**loops, "loops": {**loops["loops"]}}
    e = dict(broken["loops"]["p53"])
    e["loop_targets"] = e["loop_targets"] + ["MDM2"]
    broken["loops"]["p53"] = e
    issues = pct.validate_loops(broken)
    assert any(s == "error" and m == "p53" and "overlaps" in msg for s, m, msg in issues)


def test_v03_catalog_unchanged_by_loop_layer():
    # the loop layer must not alter the v0.3 fields used by published scores
    assert pct.load_catalog()["version"] == "0.3"
    assert "PHD2" in pct.get_genes("HIF", "feedback")


def test_get_loop_unknown_module():
    with pytest.raises(KeyError):
        pct.get_loop("NonExistentModule")
