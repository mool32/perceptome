"""9-anchor reference table — drug class × module pairs validated by Paper 4.1.

All numbers locked from PAPER4_1_CANONICAL_RESULTS.md (2026-05-09).
Three validation columns:
  block5_z, block5_q  — Block 5 v1.2 background null (1000 non-panel drugs, BH FDR)
                        PASS criterion: q < 0.10
  block2_holdout      — Block 2 independent-holdout drug list, PASSed at p<0.05
  block4_LD           — Layer Differential = |z_act| − |z_read| in Block 4 18-class panel
                        PT-PASS criterion: LD > +0.30

Roles:
  validated_rescue    — primary finding of Paper 4.1 (6 anchors)
  positive_control    — stress-axis assayed for null calibration (3 anchors); passed
                        but not claimed as primary finding

Procedure that produced these numbers, in plain words:
  For each (class, module, sign) triple, score every CMap signature on the module's
  TF-target gene panel (mean z across genes). Mean across class members = observed.
  Sample 1000 random non-panel drugs as background. Permute 10K class-sized samples
  to build a null distribution. One-sided p value vs expected sign. BH FDR across
  triples. PASS at q < 0.10.
"""

import pandas as pd

ANCHORS = pd.DataFrame([
    # ── 6 validated mechanism-pathway rescues ──────────────────────────────────
    {"class": "MEKi", "module": "ERK/MAPK", "expected_sign": -1, "role": "validated_rescue",
     "block5_z": -1.026, "block5_q": 0.0002, "block5_pass": True,
     "block2_holdout_drugs": ["PD-98059", "mek1-2-inhibitor"],
     "block2_holdout_pass": True,
     "block4_LD": 4.26, "block4_PT_pass": True,
     "notes": "Cleanest anchor; per-class scramble specificity 3× nearest class"},
    {"class": "Proteasomei", "module": "UPR-PERK", "expected_sign": +1, "role": "validated_rescue",
     "block5_z": +1.209, "block5_q": 0.0002, "block5_pass": True,
     "block2_holdout_drugs": ["ixazomib"],
     "block2_holdout_pass": True,
     "block4_LD": None, "block4_PT_pass": None,
     "notes": "Per-class scramble specificity 5× nearest class"},
    {"class": "CDKi", "module": "Cell Cycle", "expected_sign": -1, "role": "validated_rescue",
     "block5_z": None, "block5_q": None, "block5_pass": None,
     "block2_holdout_drugs": [],
     "block2_holdout_pass": None,
     "block4_LD": 3.30, "block4_PT_pass": True,
     "notes": "Block 4 layer-differential pass; not in Block 5 panel"},
    {"class": "EGFRi", "module": "ERK/MAPK", "expected_sign": -1, "role": "validated_rescue",
     "block5_z": None, "block5_q": None, "block5_pass": None,
     "block2_holdout_drugs": [],
     "block2_holdout_pass": None,
     "block4_LD": 1.07, "block4_PT_pass": True,
     "notes": "Block 4 only"},
    {"class": "PI3Ki", "module": "PI3K/PTEN", "expected_sign": +1, "role": "validated_rescue",
     "block5_z": None, "block5_q": None, "block5_pass": None,
     "block2_holdout_drugs": [],
     "block2_holdout_pass": None,
     "block4_LD": 0.81, "block4_PT_pass": True,
     "notes": "Block 4 only; FOXO targets"},
    {"class": "IKKi", "module": "NF-κB", "expected_sign": -1, "role": "validated_rescue",
     "block5_z": None, "block5_q": None, "block5_pass": None,
     "block2_holdout_drugs": [],
     "block2_holdout_pass": None,
     "block4_LD": 0.79, "block4_PT_pass": True,
     "notes": "Block 4 only"},
    # ── 3 stress-axis positive controls (Block 5 v1.2; assayed for null calibration) ──
    {"class": "HSP90i", "module": "HSF1", "expected_sign": +1, "role": "positive_control",
     "block5_z": +2.149, "block5_q": 0.0002, "block5_pass": True,
     "block2_holdout_drugs": ["NVP-AUY922", "PU-H71", "radicicol"],
     "block2_holdout_pass": True,
     "block4_LD": None, "block4_PT_pass": None,
     "notes": "Largest observed effect; clean stress-axis control"},
    {"class": "HSP90i", "module": "UPR-ATF6", "expected_sign": +1, "role": "positive_control",
     "block5_z": +0.736, "block5_q": 0.0004, "block5_pass": True,
     "block2_holdout_drugs": [],
     "block2_holdout_pass": None,
     "block4_LD": None, "block4_PT_pass": None,
     "notes": "ER-stress overlap from HSP90 inhibition"},
    {"class": "HSP90i", "module": "NRF2", "expected_sign": +1, "role": "positive_control",
     "block5_z": +0.510, "block5_q": 0.038, "block5_pass": True,
     "block2_holdout_drugs": [],
     "block2_holdout_pass": None,
     "block4_LD": None, "block4_PT_pass": None,
     "notes": "Oxidative-stress overlap"},
])

ANCHORS.attrs["paper"] = "Paper 4.1 (perceptome on drugs) closed 2026-05-09"
ANCHORS.attrs["pre_reg"] = "angle_A_prereg_v1_2_background_null.md (Block 5 v1.2)"
ANCHORS.attrs["seed"] = 42
ANCHORS.attrs["bg_pool"] = "1000 random non-panel drugs from CMap LINCS L1000 (GSE92742)"


def drug_anchors(role=None):
    """Return the 9-anchor reference table.

    Parameters
    ----------
    role : 'validated_rescue' | 'positive_control' | None
        Filter to a subset by role. None returns all 9 anchors.

    Returns
    -------
    DataFrame with columns:
        class, module, expected_sign, role,
        block5_z, block5_q, block5_pass,
        block2_holdout_drugs, block2_holdout_pass,
        block4_LD, block4_PT_pass,
        notes
    """
    if role is None:
        return ANCHORS.copy()
    return ANCHORS[ANCHORS["role"] == role].copy()
