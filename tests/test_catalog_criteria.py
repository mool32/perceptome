"""Table S9b — the four operational criteria per module, and what the paper says about them."""

import pytest

import perceptome as pct


@pytest.fixture(scope="module")
def mc():
    return pct.module_criteria()


def test_all_44_modules_are_annotated(mc):
    assert len(mc) == 44 and set(mc.index) == set(pct.list_modules())
    assert set(mc["confidence"]) <= {"rec.", "high", "medium", "low"}


def test_feedback_is_named_for_38_modules_and_unclear_for_six(mc):
    """'A feedback element is named for 38 of 44 modules; six are marked unclear'."""
    assert mc["feedback_named"].sum() == 38
    assert set(mc.index[~mc["feedback_named"]]) == {"AR", "PR", "PPARγ", "NPAS4", "PI3K/PTEN", "Hippo"}
    assert (mc.loc[~mc["feedback_named"], "feedback"].map(len) == 0).all()


def test_receptor_field_is_documented_for_28_and_blank_for_16(mc):
    """The sixteen modules whose annotation was recorded at construction (AhR ... p53)."""
    assert mc["receptor_field_recorded"].sum() == 28
    assert set(mc.index[~mc["receptor_field_recorded"]]) == {
        "AhR", "Circadian", "ERK/MAPK", "HIF", "HSF1", "Hedgehog", "Hippo", "NF-κB",
        "Notch", "SREBP", "UPR-ATF6", "VDR", "Wnt", "cGAS-STING", "mTOR", "p53"}
    assert (mc.loc[~mc["receptor_field_recorded"], "annotation"] == "recorded at catalog construction").all()


def test_named_feedback_is_absent_from_the_scored_genes_for_ten_modules(mc):
    """Note S-N / S9b: for these the score cannot be reading the feedback arm (the paper counts 11,
    including Hippo's withdrawn entry)."""
    absent = set(mc.index[mc["feedback_in_scored_genes"] == False])    # noqa: E712 (NaN-aware comparison)
    assert absent == {"AMPK", "Autophagy", "Calcium", "Cell Cycle", "ER", "Notch", "SREBP", "TR", "cGAS-STING", "mTOR"}
    assert mc["feedback_in_scored_genes"].isna().sum() == 6


def test_recorded_modules_agree_with_the_catalog_role_lists(mc):
    """The recorded annotations of Table S9b are the role lists the catalog has carried since construction."""
    cat = pct.load_catalog()["modules"]
    recorded = set(mc.index[mc["annotation"] == "recorded at catalog construction"])
    annotated = {m for m, v in cat.items() if "sensor_genes" in v}
    assert annotated - recorded == {"NPAS4"}          # re-curated in the paper (feedback unclear, low confidence)
    assert recorded <= annotated
    for m in recorded:
        for printed, key in ((mc.loc[m, "sensor"], "sensor_genes"), (mc.loc[m, "feedback"], "feedback_genes")):
            if printed:
                assert printed == cat[m][key][:len(printed)], (m, key)


def test_hippo_feedback_entry_is_withdrawn():
    """Table S9b: NKD1/NKD2 are Wnt antagonists and the BIRCs are YAP targets, not feedback inhibitors."""
    assert pct.get_genes("Hippo", "feedback") == []
    assert "feedback_note" in pct.get_module_info("Hippo")
    assert "NKD1" in pct.get_genes("Wnt", "feedback")
    assert pct.validate_catalog() == []
