"""Build the v0.4 loop layer (perceptome/catalog/data/loops_v04.json).

The loop layer sits NEXT TO the v0.3 module catalog; it does not change any
v0.3 field, so every published score (readiness, activity, eigenspace, HPA
reference) is reproduced bit-for-bit.

What it adds, per module:
  - feedback_core       2-5 negative-feedback genes that are themselves induced
                        by the module (transcriptional feedback), curated from
                        primary literature (see `refs`).
  - loop_targets        module output genes with the feedback genes REMOVED, so a
                        target-vs-feedback coupling metric never correlates a gene
                        with itself (in v0.3, 15 of 17 feedback lists overlapped
                        the activity genes).
  - feedback_mechanism  how the loop closes.
  - feedback_is_output  True when the module's output IS its feedback (chaperone
                        titration, clock, ligand catabolism): the coupling metric
                        is undefined there, only noise-compression applies.
  - feedback_readable   False when the feedback is post-translational only
                        (e.g. KEAP1 for NRF2): not measurable in mRNA.
  - perturbations       knockdown classes with opposite predictions:
                          break       remove the feedback  -> output up, cell-to-cell noise up
                          input_down  remove the input     -> output down
                          input_up    raise the input with the loop intact -> output up, noise compressed
  - tier                A  active in standard screening lines (K562 / RPE1) with clean perturbations
                        B  testable in a specific context (named in `context`)
                        C  needs an exogenous stimulus
                        D  needs a ligand / lineage absent from standard lines
                        E  no clean transcriptional loop -> not a loop-integrity candidate
  - context, v03_issues, refs

Status: CURATED, NOT YET VALIDATED. Tiers are expected from cell-line biology
and have not been checked against the Replogle et al. 2022 perturbation list.

Run:  python scripts/12_build_loops_v04.py
"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V03 = ROOT / "perceptome" / "catalog" / "data" / "modules_v03.json"
OUT = ROOT / "perceptome" / "catalog" / "data" / "loops_v04.json"

# Symbols used in v0.3 that are aliases or superseded HGNC symbols.
ALIASES = {
    "PHD2": "EGLN1",
    "GILZ": "TSC22D3",
    "CTGF": "CCN2",
    "CYR61": "CCN1",
    "G6PC": "G6PC1",
    "ACPP": "ACP3",
}


def P(brk=(), down=(), up=()):
    return {"break": list(brk), "input_down": list(down), "input_up": list(up)}


LOOPS = {
    # ---------------- tier A: active in K562 / RPE1, clean perturbations ----------------
    "ERK/MAPK": dict(
        tier="A",
        feedback_core=["DUSP6", "DUSP4", "SPRY2", "SPRY4", "SPRED1"],
        loop_targets=["FOS", "EGR1", "ETV4", "ETV5", "PHLDA1"],
        feedback_mechanism="MEK/ERK-induced phosphatases (DUSP4/6) and RTK-RAS inhibitors (SPRY/SPRED)",
        perturbations=P(brk=["DUSP6", "SPRY2", "SPRY4", "SPRED1"],
                        down=["GRB2", "SOS1", "RAF1", "MAP2K1", "MAPK1"], up=["NF1"]),
        context="K562: pathway driven by BCR-ABL",
        v03_issues="v0.3 feedback list has 21 genes incl. GADD45A/B/G (shared with NF-kB); DUSP1/DUSP6 were also targets. "
                   "DUSP1 dropped from both sets (induced by GR/p38 as much as by ERK). PHLDA1 added as MEK-dependent output.",
        refs=["Pratilas et al. 2009 PNAS", "Legewie et al. 2009 Mol Syst Biol"],
    ),
    "JAK-STAT": dict(
        tier="A",
        feedback_core=["SOCS1", "SOCS3", "CISH"],
        loop_targets=["IRF1", "BCL2L1", "MYC", "PIM1", "OSMR"],
        feedback_mechanism="STAT-induced SOCS/CIS proteins inhibit JAKs and compete for receptor docking sites",
        perturbations=P(brk=["SOCS3", "SOCS1", "CISH", "PTPN2"],
                        down=["STAT5A", "STAT5B", "STAT3", "JAK2", "JAK1"]),
        context="K562: constitutive STAT5 via BCR-ABL (CISH is the canonical STAT5 feedback); "
                "targets are STAT3/STAT1-oriented, PIM1/BCL2L1/MYC are shared STAT5 targets",
        v03_issues="no feedback_genes in v0.3; SOCS1/SOCS3 sat in activity_genes",
        refs=["Starr et al. 1997 Nature", "Matsumoto et al. 1997 Blood (CIS/STAT5)"],
    ),
    "HIF": dict(
        tier="A",
        feedback_core=["EGLN3"],
        loop_targets=["VEGFA", "SLC2A1", "LDHA", "PGK1", "BNIP3", "CA9"],
        feedback_mechanism="HIF induces the prolyl hydroxylase EGLN3 (PHD3); EGLN1/VHL degradation machinery is constitutive",
        perturbations=P(brk=["VHL", "EGLN1"], down=["HIF1A", "ARNT"]),
        context="normoxic lines: basal activity low, VHL/EGLN1 knockdown strongly induces",
        v03_issues="'PHD2' is an alias of EGLN1 (duplicate entry); EGLN3 was also an activity gene",
        refs=["Epstein et al. 2001 Cell", "Aprelikova et al. 2004 J Cell Biochem"],
    ),
    "HSF1": dict(
        tier="A",
        feedback_core=["HSPA1A", "HSPA1B", "DNAJB1", "HSP90AA1"],
        loop_targets=["HSPA6", "HSPB1", "BAG3"],
        feedback_mechanism="chaperone titration: HSF1-induced HSP70/HSP90 rebind and inactivate HSF1",
        feedback_is_output=True,
        perturbations=P(brk=["HSP90AB1", "HSP90AA1", "HSPA8"], down=["HSF1"]),
        context="chaperone knockdowns are also proteotoxic; compare against HSF1 knockdown",
        v03_issues="feedback and activity lists overlapped (HSPA1A, HSP90AA1); module internally incoherent "
                   "in the step10 audit (gene coherence -0.06, FRAGILE)",
        refs=["Zou et al. 1998 Cell", "Zheng et al. 2016 eLife", "Krakowiak et al. 2018 eLife"],
    ),
    "UPR-ATF6": dict(
        tier="A",
        feedback_core=["HSPA5"],
        loop_targets=["CALR", "PDIA4", "PDIA6", "HSP90B1"],
        feedback_mechanism="BiP (HSPA5) is both the ER-stress sensor and an ATF6-induced output that re-represses ATF6",
        feedback_is_output=True,
        perturbations=P(brk=["HSPA5"], down=["ATF6", "MBTPS1", "MBTPS2"]),
        context="HSPA5 knockdown activates all three UPR arms; ATF6 knockdown gives arm specificity; "
                "MBTPS1/2 are shared with SREBP (built-in specificity test)",
        v03_issues="HSPA5 and CALR were both feedback and activity genes; HSPA5/PDIA4 shared with UPR-IRE1 in paper-4 sets",
        refs=["Bertolotti et al. 2000 Nat Cell Biol", "Shen et al. 2002 Dev Cell"],
    ),
    "UPR-PERK": dict(
        tier="A",
        feedback_core=["PPP1R15A"],
        loop_targets=["DDIT3", "ASNS", "ATF3", "TRIB3", "SESN2"],
        feedback_mechanism="ATF4/CHOP induce GADD34 (PPP1R15A), which dephosphorylates eIF2a",
        perturbations=P(brk=["PPP1R15A"], down=["EIF2AK3", "ATF4"]),
        context="readout is the integrated stress response; in erythroid K562 the eIF2a kinase HRI (EIF2AK1) "
                "is a second input, so EIF2AK3 knockdown alone may under-shoot",
        v03_issues="no feedback_genes in v0.3; PPP1R15A sat in activity_genes",
        refs=["Novoa et al. 2001 J Cell Biol", "Marciniak et al. 2004 Genes Dev"],
    ),
    "NRF2": dict(
        tier="A",
        feedback_core=[],
        loop_targets=["NQO1", "HMOX1", "TXNRD1", "GCLC", "GCLM", "SLC7A11", "GPX2", "FTH1"],
        feedback_mechanism="post-translational: KEAP1-CUL3 ubiquitinates NRF2; no reliable transcriptional feedback",
        feedback_readable=False,
        perturbations=P(brk=["KEAP1", "CUL3"], down=["NFE2L2"]),
        context="break-only test: coupling metric not applicable, noise-compression metric applies",
        v03_issues="no feedback_genes in v0.3",
        refs=["Itoh et al. 1999 Genes Dev", "Kobayashi et al. 2004 Mol Cell Biol"],
    ),
    "SREBP": dict(
        tier="A",
        feedback_core=["INSIG1"],
        loop_targets=["HMGCR", "FASN", "SCD", "ACLY", "LDLR"],
        feedback_mechanism="SREBP induces INSIG1, which retains SCAP-SREBP in the ER",
        perturbations=P(brk=["INSIG1", "INSIG2"], down=["SCAP", "SREBF2", "SREBF1", "MBTPS1", "MBTPS2"]),
        context="MBTPS1/2 shared with UPR-ATF6 (specificity test); targets mix SREBP1 lipogenic and SREBP2 sterol arms",
        v03_issues="v0.3 feedback included LDLR/PCSK9 (uptake arm, metabolite-mediated) and LDLR was also a target",
        refs=["Yang et al. 2002 Cell", "Horton et al. 2002 J Clin Invest"],
    ),
    # ---------------- tier B: testable in a specific context ----------------
    "p53": dict(
        tier="B",
        feedback_core=["MDM2", "PPM1D"],
        loop_targets=["CDKN1A", "BAX", "BBC3", "GADD45A", "PMAIP1", "SESN2", "TIGAR", "DDB2", "RRM2B"],
        feedback_mechanism="p53 induces its E3 ligase MDM2 and the phosphatase WIP1 (PPM1D)",
        perturbations=P(brk=["MDM2", "PPM1D", "MDM4"], down=["TP53", "ATM"]),
        context="requires TP53-wild-type cells (RPE1); K562 is TP53-null",
        v03_issues="MDM2 was both feedback and activity gene; PPP1R13L shared with NF-kB",
        refs=["Lahav et al. 2004 Nat Genet", "Batchelor et al. 2008 Mol Cell"],
    ),
    "UPR-IRE1": dict(
        tier="B",
        feedback_core=["DNAJB9"],
        loop_targets=["SEC61A1", "EDEM1", "ERO1A"],
        feedback_mechanism="XBP1s induces ERdj4 (DNAJB9), which recruits BiP to repress IRE1",
        perturbations=P(brk=["DNAJB9", "HSPA5"], down=["ERN1", "XBP1"]),
        context="only three unshared targets (HSPA5 removed: shared with ATF6)",
        v03_issues="no feedback_genes in v0.3; DNAJB9 sat in activity_genes",
        refs=["Amin-Wetzel et al. 2017 Cell"],
    ),
    "mTOR": dict(
        tier="B",
        feedback_core=[],
        loop_targets=["SLC7A5", "SREBF1", "VEGFA", "HIF1A", "RPS6", "EIF4EBP1"],
        feedback_mechanism="post-translational (S6K1-IRS1, GRB10, DEPTOR); mTORC1 output is mostly translational, "
                           "so mRNA targets are a weak readout",
        feedback_readable=False,
        perturbations=P(down=["RPTOR", "MTOR", "RHEB"], up=["TSC1", "TSC2"]),
        context="input_up (TSC1/2) vs input_down only; no clean break perturbation",
        v03_issues="v0.3 feedback list (DDIT4, SESN2, ULK1/2, DEPTOR, GRB10) mixed upstream inhibitors with mTOR substrates",
        refs=["Harrington et al. 2004 J Cell Biol", "Hsu et al. 2011 Science"],
    ),
    "Hippo": dict(
        tier="B",
        feedback_core=["LATS2", "NF2", "AMOTL2"],
        loop_targets=["CCN2", "CCN1", "ANKRD1"],
        feedback_mechanism="YAP/TAZ-TEAD induce LATS2, NF2 and AMOTL2, which re-inactivate YAP/TAZ",
        perturbations=P(brk=["LATS2", "NF2", "AMOTL2"], down=["YAP1", "WWTR1", "TEAD1"], up=["LATS1"]),
        context="adherent cells (RPE1); YAP activity low in suspension K562",
        v03_issues="v0.3 feedback held NKD1/NKD2 (Wnt) and BIRC2/3 (NF-kB); BIRC5 dropped from targets "
                   "(survivin tracks proliferation); CTGF/CYR61 renamed CCN2/CCN1",
        refs=["Moroishi et al. 2015 Genes Dev"],
    ),
    # ---------------- tier C: needs an exogenous stimulus ----------------
    "NF-κB": dict(
        tier="C",
        feedback_core=["NFKBIA", "TNFAIP3", "NFKBIE"],
        loop_targets=["IL6", "TNF", "CXCL10", "CCL2", "ICAM1", "IL1B", "CXCL8", "BCL2L1", "MMP9", "PTGS2"],
        feedback_mechanism="NF-kB resynthesises its inhibitors IkBa/IkBe and the deubiquitinase A20",
        perturbations=P(brk=["NFKBIA", "TNFAIP3"], down=["RELA", "IKBKB", "IKBKG"]),
        context="resting K562/RPE1 have little NF-kB activity; use TNF/IL-1 stimulation data",
        v03_issues="v0.3 feedback list had 15 genes incl. UBE2I, PIAS4, GADD45A/B/G, BIRC2/3; NFKBIA/TNFAIP3 were also targets",
        refs=["Hoffmann et al. 2002 Science", "Lee et al. 2000 Science", "Nelson et al. 2004 Science"],
    ),
    "Wnt": dict(
        tier="C",
        feedback_core=["AXIN2", "NKD1", "RNF43", "ZNRF3"],
        loop_targets=["LEF1", "TCF7", "MYC", "CCND1", "SP5"],
        feedback_mechanism="beta-catenin induces AXIN2/NKD1 (destruction complex) and RNF43/ZNRF3 (receptor clearance)",
        perturbations=P(brk=["AXIN2", "RNF43", "ZNRF3"], down=["CTNNB1", "TCF7L2"], up=["APC", "CSNK1A1", "GSK3B"]),
        context="low basal Wnt in K562/RPE1; MYC/CCND1 also track proliferation",
        v03_issues="v0.3 feedback had 23 genes mixing extracellular antagonists (DKK/SFRP/WIF1) with cell-intrinsic feedback",
        refs=["Lustig et al. 2002 Mol Cell Biol", "Koo et al. 2012 Nature", "Hao et al. 2012 Nature"],
    ),
    "TGF-β": dict(
        tier="C",
        feedback_core=["SMAD7", "SKIL", "PMEPA1"],
        loop_targets=["SERPINE1", "COL1A1", "CCN2", "FN1", "TGFBI"],
        feedback_mechanism="SMAD-induced inhibitory SMAD7, SnoN (SKIL) and PMEPA1",
        perturbations=P(brk=["SMAD7", "SKIL"], down=["SMAD4", "SMAD3", "TGFBR1", "TGFBR2"]),
        context="needs TGF-b stimulation or autocrine signalling",
        v03_issues="no feedback_genes in v0.3; SMAD7 sat in activity_genes; CTGF renamed CCN2",
        refs=["Nakao et al. 1997 Nature", "Stroschein et al. 1999 Science"],
    ),
    "BMP": dict(
        tier="C",
        feedback_core=["SMAD6", "BAMBI"],
        loop_targets=["ID1", "ID2", "ID3"],
        feedback_mechanism="BMP-induced inhibitory SMAD6 and pseudoreceptor BAMBI",
        perturbations=P(brk=["SMAD6"], down=["SMAD1", "SMAD5", "BMPR1A", "SMAD4"]),
        context="needs BMP stimulation",
        v03_issues="no feedback_genes in v0.3; SMAD6/BAMBI sat in activity_genes",
        refs=["Imamura et al. 1997 Nature", "Onichtchouk et al. 1999 Nature"],
    ),
    "Type I IFN": dict(
        tier="C",
        feedback_core=["USP18"],
        loop_targets=["ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IFI44L"],
        feedback_mechanism="ISG USP18 displaces JAK1 from IFNAR2",
        perturbations=P(brk=["USP18"], down=["IFNAR1", "IFNAR2", "STAT1", "STAT2", "IRF9"]),
        context="needs IFN or tonic IFN signalling",
        v03_issues="no feedback_genes in v0.3; USP18 sat in activity_genes",
        refs=["Malakhova et al. 2006 EMBO J", "Sarasin-Filipowicz et al. 2009 Mol Cell Biol"],
    ),
    "cGAS-STING": dict(
        tier="C",
        feedback_core=["TREX1"],
        loop_targets=["IFNB1", "CXCL10", "ISG15", "IFIT1", "MX1", "OAS1"],
        feedback_mechanism="IRF3/IFN-inducible exonuclease TREX1 degrades the cytosolic DNA trigger",
        perturbations=P(brk=["TREX1"], down=["CGAS", "STING1", "TBK1", "IRF3"]),
        context="readout overlaps Type I IFN; specificity only through input knockdowns",
        v03_issues="v0.3 feedback mixed NFKBIA/NFKBIB (NF-kB arm) and ADAR (MDA5 arm)",
        refs=["Stetson et al. 2008 Cell", "Ablasser et al. 2014 J Immunol"],
    ),
    "cAMP/CREB": dict(
        tier="C",
        feedback_core=["CREM", "PDE4B", "PDE4D"],
        loop_targets=["BDNF", "NR4A1", "EGR1", "ARC", "FOS", "FOSB", "ATF3", "NR4A2", "NR4A3"],
        feedback_mechanism="CREB induces the repressor ICER (CREM) and cAMP-degrading PDE4",
        perturbations=P(brk=["PDE4B", "PDE4D", "CREM"], down=["CREB1", "PRKACA"]),
        context="needs cAMP stimulation; FOS/EGR1 shared with ERK, ATF3 with UPR-PERK",
        v03_issues="no feedback_genes in v0.3",
        refs=["Molina et al. 1993 Cell", "D'Sa et al. 2002 J Neurochem"],
    ),
    "AhR": dict(
        tier="C",
        feedback_core=["AHRR", "CYP1A1", "TIPARP"],
        loop_targets=["CYP1B1", "NQO1"],
        feedback_mechanism="AHR induces its repressor AHRR, the ligand-degrading CYP1A1 and the ADP-ribosyltransferase TIPARP",
        feedback_is_output=True,
        perturbations=P(brk=["AHRR", "TIPARP"], down=["AHR", "ARNT"]),
        context="needs AhR ligand; ARNT shared with HIF (specificity test); NQO1 shared with NRF2",
        v03_issues="all three v0.3 feedback genes were also activity genes",
        refs=["Mimura et al. 1999 Genes Dev", "MacPherson et al. 2013 Nucleic Acids Res"],
    ),
    # ---------------- tier D: needs a ligand / lineage ----------------
    "NFAT": dict(
        tier="D",
        feedback_core=["RCAN1"],
        loop_targets=["IL2", "CSF2", "IL4", "FASLG", "GRIA1", "KCNA4"],
        feedback_mechanism="NFAT induces RCAN1 (DSCR1), an inhibitor of calcineurin",
        perturbations=P(brk=["RCAN1"], down=["PPP3CA", "PPP3CB", "PPP3R1", "NFATC1", "NFATC2"]),
        context="needs Ca2+ stimulation in T cells / neurons",
        v03_issues="no feedback_genes in v0.3; RCAN1 sat in activity_genes (also in Calcium)",
        refs=["Yang et al. 2000 Circ Res", "Rothermel et al. 2003 Trends Cardiovasc Med"],
    ),
    "GR": dict(
        tier="D",
        feedback_core=["FKBP5"],
        loop_targets=["TSC22D3", "SGK1", "PER1", "KLF13"],
        feedback_mechanism="GR induces FKBP5, which lowers GR ligand affinity and nuclear translocation",
        perturbations=P(brk=["FKBP5"], down=["NR3C1"]),
        context="needs glucocorticoid",
        v03_issues="no feedback_genes in v0.3; FKBP5 sat in activity_genes; DUSP1 dropped (shared with ERK)",
        refs=["Vermeer et al. 2003 J Clin Endocrinol Metab", "Binder 2009 Psychoneuroendocrinology"],
    ),
    "VDR": dict(
        tier="D",
        feedback_core=["CYP24A1"],
        loop_targets=["CAMP", "TRPV6", "S100A8", "DEFB4A"],
        feedback_mechanism="ligand catabolism: VDR induces CYP24A1, which degrades 1,25(OH)2D3",
        feedback_is_output=True,
        perturbations=P(brk=["CYP24A1"], down=["VDR"]),
        context="needs vitamin D",
        v03_issues="CYP24A1 was both feedback and activity gene",
        refs=["Jones et al. 2012 Arch Biochem Biophys"],
    ),
    "RAR": dict(
        tier="D",
        feedback_core=["CYP26A1"],
        loop_targets=["CRABP2", "TGM2", "RARB", "HOXA1"],
        feedback_mechanism="ligand catabolism: RAR induces CYP26A1, which degrades retinoic acid",
        feedback_is_output=True,
        perturbations=P(brk=["CYP26A1"], down=["RARA"]),
        context="needs retinoic acid",
        v03_issues="no feedback_genes in v0.3; CYP26A1 sat in activity_genes",
        refs=["Ray et al. 1997 J Biol Chem"],
    ),
    "FXR": dict(
        tier="D",
        feedback_core=["NR0B2"],
        loop_targets=["ABCB11", "SLC51A", "FGF19", "BAAT"],
        feedback_mechanism="FXR induces the repressor SHP (NR0B2)",
        perturbations=P(brk=["NR0B2"], down=["NR1H4"]),
        context="hepatocyte / enterocyte only",
        v03_issues="no feedback_genes in v0.3; NR0B2 sat in activity_genes",
        refs=["Goodwin et al. 2000 Mol Cell"],
    ),
    "Circadian": dict(
        tier="D",
        feedback_core=["PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2"],
        loop_targets=["DBP", "TEF"],
        feedback_mechanism="the transcription-translation loop is the module itself",
        feedback_is_output=True,
        perturbations=P(brk=["CRY1", "CRY2", "PER2"], down=["CLOCK", "BMAL1"]),
        context="unsynchronised cultured cells average the clock away",
        v03_issues="CRY1/NR1D1/PER1/PER2 were both feedback and activity genes",
        refs=["Takahashi 2017 Nat Rev Genet"],
    ),
    "Hedgehog": dict(
        tier="D",
        feedback_core=["PTCH1", "HHIP"],
        loop_targets=["GLI1", "PTCH2", "FOXF1"],
        feedback_mechanism="GLI induces the receptor PTCH1 and the ligand trap HHIP",
        feedback_is_output=True,
        perturbations=P(down=["SMO", "GLI1", "GLI2"], up=["PTCH1", "SUFU"]),
        context="inactive outside developmental / ciliated contexts",
        v03_issues="all three v0.3 feedback genes were also activity genes",
        refs=["Chuang & McMahon 1999 Nature"],
    ),
    "Notch": dict(
        tier="D",
        feedback_core=["NRARP"],
        loop_targets=["HES1", "HEY1", "HEY2"],
        feedback_mechanism="Notch induces NRARP, which destabilises the NICD-RBPJ complex",
        perturbations=P(brk=["NRARP"], down=["RBPJ", "NOTCH1", "MAML1"]),
        context="needs ligand-expressing neighbours",
        v03_issues="NRARP/DTX1 were both feedback and activity genes",
        refs=["Lamar et al. 2001 Genes Dev"],
    ),
}

# Modules with no clean transcriptional negative-feedback loop.
TIER_E = {
    "AMPK": "energy sensing is post-translational; no transcriptional feedback curated",
    "Autophagy": "flux is post-translational; transcriptional TFEB loop not curated",
    "PI3K/PTEN": "feedback is post-translational (S6K-IRS1, AKT-FOXO)",
    "Insulin/FOXO": "feedback is post-translational",
    "Cell Cycle": "oscillator, not a homeostatic negative-feedback module",
    "NPAS4": "post-translational turnover (see v0.3 note)",
    "Calcium": "loop curated under NFAT (RCAN1)",
    "AR": "no ligand in standard lines; FKBP5 feedback curated under GR",
    "MR": "no ligand in standard lines",
    "PR": "no ligand in standard lines",
    "ER": "no ligand in standard lines",
    "PXR/CAR": "no ligand in standard lines",
    "LXR": "no clean transcriptional feedback on the receptor",
    "PPARα": "no clean transcriptional feedback on the receptor",
    "PPARγ": "no clean transcriptional feedback on the receptor",
    "TR": "no clean transcriptional feedback on the receptor",
}


def main():
    v03 = json.loads(V03.read_text())["modules"]
    loops = {}
    for name, d in LOOPS.items():
        assert name in v03, f"{name} not in v0.3 catalog"
        entry = {
            "tier": d["tier"],
            "feedback_core": d["feedback_core"],
            "loop_targets": d["loop_targets"],
            "feedback_mechanism": d["feedback_mechanism"],
            "feedback_is_output": d.get("feedback_is_output", False),
            "feedback_readable": d.get("feedback_readable", True),
            "perturbations": d["perturbations"],
            "context": d["context"],
            "v03_issues": d["v03_issues"],
            "refs": d["refs"],
        }
        loops[name] = entry
    for name, why in TIER_E.items():
        assert name in v03, f"{name} not in v0.3 catalog"
        loops[name] = {
            "tier": "E", "feedback_core": [], "loop_targets": [],
            "feedback_mechanism": why, "feedback_is_output": False, "feedback_readable": False,
            "perturbations": P(), "context": "", "v03_issues": "", "refs": [],
        }
    missing = set(v03) - set(loops)
    assert not missing, f"modules without a loop entry: {missing}"

    out = {
        "version": "0.4-loops",
        "base_catalog": "0.3",
        "status": "curated, not yet validated; tiers expected from cell-line biology, "
                  "not yet checked against the Replogle et al. 2022 perturbation list",
        "n_modules": len(loops),
        "aliases": ALIASES,
        "tiers": {
            "A": "active in standard screening lines (K562 / RPE1) with clean perturbations",
            "B": "testable in a specific context (see context)",
            "C": "needs an exogenous stimulus",
            "D": "needs a ligand / lineage absent from standard lines",
            "E": "no clean transcriptional loop; not a loop-integrity candidate",
        },
        "perturbation_classes": {
            "break": "remove the feedback -> output up, cell-to-cell noise up, target-feedback coupling lost",
            "input_down": "remove the input -> output down",
            "input_up": "raise the input with the loop intact -> output up, noise compressed",
        },
        "loops": dict(sorted(loops.items())),
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    from collections import Counter
    print(f"Wrote {OUT.relative_to(ROOT)}: {len(loops)} modules, tiers {dict(sorted(Counter(l['tier'] for l in loops.values()).items()))}")


if __name__ == "__main__":
    main()
