"""Extract the Supporting-Information tables the tool bundles from the paper's text layer.

Source: Spiro (2026), "Diverse cancers converge on one region of normal cell-state space in a
cancer-blind coordinate system", Supporting Information (Zenodo 10.5281/zenodo.20542130).

Input is the text layer of the PDF:

    pdftotext -layout paper3_v11_FULL.pdf paper_layout.txt
    python scripts/14_extract_paper_tables.py paper_layout.txt

Output (perceptome/convergence/data/):
    tableS1_datasets.csv       the 25 cancer cohorts
    tableS2_ranking154.csv     cross-cancer recurrence / depth / occupancy for the 154 normal cell types
    tableS4_modules.csv        per-module cancer -> core contributions (44 modules, 25 cohorts)
    tableS6_complexity.csv     transcriptome complexity of the cells assigned to each normal cell type
    tableS7_core_signature.csv the core's target-defining module signature
    tableS8_origins.csv        cancer -> cell of origin
and, into perceptome/catalog/data/:
    module_criteria_s9b.json   Table S9b, the four operational criteria per module (as printed)

Each table is verified against a total the paper states itself (448,399 QC cells; 574,451 assignments;
contributions summing to 0.711; 154 / 44 / 25 rows) before it is written; the script exits non-zero
if a check fails, so a transcription slip cannot go unnoticed.
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "perceptome" / "convergence" / "data"
CATALOG_DATA = ROOT / "perceptome" / "catalog" / "data"

_MATH = str.maketrans({"\U0001D705": "κ", "\U0001D6FC": "α", "\U0001D6FE": "γ", "\U0001D6FD": "β"})
_NUM = r"-?\d+(?:\.\d+)?"


def norm(s):
    """Math-italic letters used in the PDF -> the Greek letters used by the catalog."""
    return s.translate(_MATH).strip()


def section(lines, start_pat, end_pat):
    """Lines from the first match of start_pat up to (excluding) the next match of end_pat."""
    i = next(k for k, l in enumerate(lines) if re.match(start_pat, l))
    j = next((k for k in range(i + 1, len(lines)) if re.match(end_pat, lines[k])), len(lines))
    return lines[i:j]


# ── Table S1: hand-entered (cells wrap over several lines in the PDF) ─────────

S1 = [
    ("HCC", "Sun 2021", "liver", "10x", 54978, "epithelial", "author"),
    ("BRCA", "Wu 2021", "breast", "10x", 65104, "epithelial", "author"),
    ("Glioblastoma", "Census", "brain", "10x", 25661, "glial", "Census"),
    ("Ovarian serous", "Census", "ovary", "10x", 28950, "epithelial", "Census"),
    ("RCC", "Census", "kidney", "10x", 26012, "epithelial", "Census"),
    ("Melanoma", "Census", "skin", "10x", 12030, "neural-crest", "Census"),
    ("Prostate", "Census", "prostate", "10x", 26109, "epithelial", "Census"),
    ("HNSCC", "Census", "head/neck", "10x", 14808, "epithelial", "author"),
    ("Neuroblastoma", "Census", "adrenal", "10x", 30000, "neural-crest", "Census"),
    ("B-ALL", "Census", "bone marrow", "10x", 29132, "hematopoietic", "Census"),
    ("DLBCL", "Census", "lymph node", "10x", 12786, "hematopoietic", "Census"),
    ("Squamous lung", "Census", "lung", "10x", 19707, "epithelial", "Census"),
    ("TNBC", "Census", "breast", "10x", 15419, "epithelial", "Census"),
    ("PDAC", "Peng 2019", "pancreas", "10x", 505, "epithelial", "author"),
    ("NSCLC", "Census", "lung", "10x", 204, "epithelial", "Census"),
    ("AML", "van Galen 2019", "bone marrow", "Seq-Well", 6274, "hematopoietic", "author"),
    ("MM", "Ledergor 2018", "bone marrow", "InDrops", 5138, "hematopoietic", "author"),
    ("BCC", "Yost 2019", "skin", "10x", 5526, "epithelial", "heuristic"),
    ("Colorectal", "Pelka 2021", "colon", "10x/Seq-Well", 23712, "epithelial", "author"),
    ("Small-cell lung", "Census", "lung", "10x", 22388, "neuroendocrine", "Census"),
    ("Gastric", "Census", "stomach", "10x", 4632, "epithelial", "Census"),
    ("Germ-cell (embryonal)", "GSE228501", "testis", "10x", 443, "germ cell", "marker"),
    ("Osteosarcoma", "GSE162454", "bone", "10x", 16214, "mesenchymal", "marker"),
    ("Ewing sarcoma", "GSE130024", "bone (PDX)", "10x", 2007, "mesenchymal", "marker"),
    ("Seminoma", "GSE197778", "testis", "10x", 660, "germ cell", "marker"),
]


def table_s1(lines):
    df = pd.DataFrame(S1, columns=["cancer", "source", "tissue", "platform", "n_cells_QC",
                                   "origin_class", "malignant_label"])
    # every cancer name and cell count must also be present in the PDF text
    text = "\n".join(section(lines, r"^Table S1\. Datasets", r"^Table S2\."))
    for r in df.itertuples():
        assert str(r.n_cells_QC) in text, f"S1: {r.cancer} {r.n_cells_QC} not found in the PDF text"
    assert len(df) == 25
    assert df["n_cells_QC"].sum() == 448399, df["n_cells_QC"].sum()   # 'full quality-controlled pool ... 448,399'
    return df


# ── Table S2: 154-cell ranking ───────────────────────────────────────────────

_S2_ROW = re.compile(rf"^(\S.*?)\s{{2,}}({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+candidate\s+({_NUM})\s*$")
_S2_QUAR = re.compile(rf"^(erythrocytes|platelets)\s{{2,}}quarantined\s+({_NUM})\s*$")


def table_s2(lines):
    block = section(lines, r"^Table S2\. Full 154-cell ranking", r"^\s*• Table S3")
    rows, prev = [], None
    for l in block:
        m, q = _S2_ROW.match(l), _S2_QUAR.match(l)
        if m:
            rows.append({"cell_type": m.group(1).strip(), "recurrence": int(m.group(2)),
                         "median_depth_pct": float(m.group(3)), "mean_depth_pct": float(m.group(4)),
                         "occupancy_pct": float(m.group(5)), "qc_status": "candidate",
                         "pre_quarantine_occupancy_pct": float(m.group(6))})
            prev = rows[-1]
        elif q:
            rows.append({"cell_type": q.group(1), "recurrence": None, "median_depth_pct": None,
                         "mean_depth_pct": None, "occupancy_pct": None, "qc_status": "quarantined",
                         "pre_quarantine_occupancy_pct": float(q.group(2))})
            prev = None
        elif prev is not None and re.match(r"^[a-z][a-z -]*$", l.strip()) and len(l) - len(l.lstrip()) < 4 and l.strip():
            prev["cell_type"] += " " + l.strip()   # wrapped remainder of a long name ('... absorptive' / 'cells')
            prev = None
    df = pd.DataFrame(rows)
    assert len(df) == 154, len(df)
    assert df["cell_type"].is_unique
    cand = df[df.qc_status == "candidate"]
    assert len(cand) == 152
    # occupancy over the 152 candidates is a share of all malignant cells
    assert abs(cand["occupancy_pct"].sum() - 100.0) < 0.5, cand["occupancy_pct"].sum()
    top = cand.sort_values("recurrence", ascending=False).iloc[0]
    assert top["recurrence"] == 18
    return df


# ── Table S4: per-module cancer -> core contributions ────────────────────────

_S4_ROW = re.compile(rf"^(\S.*?)\s{{2,}}({_NUM})\s+({_NUM})\s+({_NUM})\s+(25)\s*$")


def table_s4(lines):
    block = section(lines, r"^\s*• Table S3 \(inline", r"^\s*• Table S5")
    rows = []
    for l in block:
        m = _S4_ROW.match(l)
        if m and not l.startswith("module"):
            rows.append({"module": norm(m.group(1)), "mean_contribution": float(m.group(2)),
                         "sd": float(m.group(3)), "pct_positive": float(m.group(4)), "n_cohorts": int(m.group(5))})
    df = pd.DataFrame(rows)
    assert len(df) == 44, len(df)
    # 'the 44 contributions sum to that cosine — here 0.711'
    assert abs(df["mean_contribution"].sum() - 0.711) < 0.002, df["mean_contribution"].sum()
    return df


# ── Table S6: transcriptome complexity per assigned cell type ────────────────

_S6_ROW = re.compile(rf"^(\S.*?)\s{{2,}}(\d+)\s+({_NUM}|nan)\s+({_NUM}|nan)\s+(candidate|quarantined)(.*)$")


def table_s6(lines):
    block = section(lines, r"^\s*\u2022 Table S5 \(inline", r"^Table S7\.")
    rows, prev = [], None
    for l in block:
        m = _S6_ROW.match(l)
        if m:
            tail = m.group(6).strip()
            rows.append({"cell_type": m.group(1).strip(), "n_cells_assigned": int(m.group(2)),
                         "median_ngenes": None if m.group(3) == "nan" else float(m.group(3)),
                         "mean_ngenes": None if m.group(4) == "nan" else float(m.group(4)),
                         "qc_status": m.group(5), "note": "" if tail in ("", "nan") else tail})
            prev = rows[-1]
        elif prev is not None and l.strip() and not re.fullmatch(r"\d+", l.strip()) \
                and not l.lstrip().startswith(("n_cells", "cell ", "Table")):   # (digits alone = page number)
            # continuation of a wrapped cell: name remainder at the left margin, note remainder far right
            pieces = re.split(r"\s{2,}", l.strip())
            if len(l) - len(l.lstrip()) < 10:
                prev["cell_type"] += " " + pieces[0]
                pieces = pieces[1:]
            if pieces and prev["qc_status"] == "candidate":
                prev["note"] = (prev["note"] + " " + " ".join(pieces)).strip()
            prev = None
    for r in rows:
        if r["qc_status"] == "quarantined":
            r["note"] = "post-cellular, low-complexity"
    df = pd.DataFrame(rows)
    assert len(df) == 154, len(df)
    assert df["cell_type"].is_unique
    assert df["n_cells_assigned"].sum() == 574451, df["n_cells_assigned"].sum()   # 'Table S6 tabulates 574,451 assignments'
    assert (df["note"] == "no cells assigned").sum() == df["median_ngenes"].isna().sum() == 10
    return df


# ── Table S7: the core's target-defining signature ───────────────────────────

_S7_ROW = re.compile(rf"^(\S.*?)\s{{2,}}({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s*$")


def table_s7(lines):
    block = section(lines, r"^Table S7\.", r"^Table S8\.")
    rows = []
    for l in block:
        m = _S7_ROW.match(l)
        if m:
            rows.append({"module": norm(m.group(1)), "target_distinctiveness": float(m.group(2)),
                         "core_suppression": float(m.group(3)), "megakaryocyte_z": float(m.group(4)),
                         "cytotrophoblast_z": float(m.group(5)), "cancer_contribution": float(m.group(6))})
    df = pd.DataFrame(rows)
    assert len(df) == 44, len(df)
    assert (df["core_suppression"] + df["target_distinctiveness"]).abs().max() < 1e-9   # suppression = -distinctiveness
    return df


# ── Table S8: cancer -> cell of origin ───────────────────────────────────────

def table_s8(lines):
    block = section(lines, r"^Table S8\.", r"^Table S9\.")
    rows = []
    for l in block:
        m = re.match(r"^\s{3,}(\S.*?)\s{2,}(\S.*?)\s{2,}(\S.*?)\s*$", l)
        if m and m.group(1) != "cancer":
            rows.append({"cancer": m.group(1), "cell_of_origin": m.group(2), "origin_class": m.group(3)})
    df = pd.DataFrame(rows)
    assert len(df) == 25, len(df)
    return df


# ── Table S9b: the four operational criteria per module ──────────────────────

def _columns(header):
    """Start offsets of the six columns of a Table S9b page, from its header line."""
    header = header.replace("\x0c", "")          # the page-break form feed precedes the header of page 2
    return [header.index(k) for k in ("module", "sensor", "TF output", "feedback element", "receptor field", "conf.")]


def _split_cells(line, cols):
    """Slice a table line at the column offsets, nudging a boundary left when it would cut a token.

    The text layer is not a perfect grid: a long cell can run a character or two into the next
    column (e.g. 'NFE2L2 heterodimerized KEAP1'), and a dash in the receptor column can sit just
    left of its header. Both are repaired here.
    """
    bounds = list(cols[1:]) + [len(line) + 1]
    adj = []
    for b in bounds:
        while 0 < b < len(line) and line[b - 1] != " " and line[b] != " ":
            b -= 1
        adj.append(b)
    cells = [line[a:b].strip() for a, b in zip([cols[0]] + adj[:-1], adj)]
    for k in (1, 2, 3):                         # a dash glued to the end of the cell before the receptor field
        if cells[k].endswith("\u2014") and len(cells[k]) > 1 and k == 3:
            cells[3] = cells[3][:-1].strip()
            cells[4] = "\u2014"
    return cells


def table_s9b(lines):
    block = section(lines, r"^Table S9b\.", r"^Note S-N\.")
    heads = [i for i, l in enumerate(block) if re.match(r"^\s?module\s+sensor\s+TF output", l)]
    assert len(heads) == 2, heads
    rows = {}
    for h, end in zip(heads, heads[1:] + [len(block)]):
        cols = _columns(block[h])
        cur = None
        for l in block[h + 1:end]:
            if not l.strip() or re.fullmatch(r"\s*\d+\s*", l):
                continue
            cells = _split_cells(l, cols)
            if cells[0]:
                cur = norm(cells[0])
                rows[cur] = [[c] if c else [] for c in cells[1:]]
            elif cur:
                for k, c in enumerate(cells[1:]):
                    if c:
                        rows[cur][k].append(c)
    out = {}
    for m, parts in rows.items():
        sensor, tf, fb, rf, conf = (" ".join(p) for p in parts)
        def genes(t):
            """Comma-separated symbols; a trailing fragment cut off by the print ('TRIB3 p\u2026') is dropped."""
            items = [re.sub(r"\s*\S*\u2026\s*$", "", g).strip() for g in t.split(",")]
            return [g for g in items if g]
        out[m] = {
            "sensor": genes(sensor),
            "tf_output": tf,
            "tf_output_truncated": tf.endswith("\u2026"),
            "feedback": [] if fb == "unclear" else genes(fb),
            "feedback_named": fb != "unclear",
            "receptor_field_recorded": rf not in ("\u2014", ""),
            "confidence": conf,
        }
    assert len(out) == 44, len(out)
    # statements of the paper itself
    assert sum(v["feedback_named"] for v in out.values()) == 38                       # 'named for 38 of 44 modules'
    assert {m for m, v in out.items() if not v["feedback_named"]} == \
        {"AR", "PR", "PPAR\u03b3", "NPAS4", "PI3K/PTEN", "Hippo"}                        # 'six are marked unclear'
    assert sum(not v["receptor_field_recorded"] for v in out.values()) == 16          # 'blank for sixteen modules'
    cat = json.loads((CATALOG_DATA / "modules_v03.json").read_text(encoding="utf-8"))["modules"]
    absent = {m for m, v in out.items() if v["feedback"] and not set(v["feedback"]) & set(cat[m]["core_genes"])}
    assert absent == {"AMPK", "Autophagy", "Calcium", "Cell Cycle", "ER", "Notch", "SREBP", "TR", "cGAS-STING", "mTOR"}, absent
    for v in out.values():       # 16 recorded at construction (blank receptor field), 28 re-curated
        v["annotation"] = "recorded at catalog construction" if not v["receptor_field_recorded"] else "re-curated"
    return out


def main(path):
    lines = Path(path).read_text(encoding="utf-8").split("\n")
    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "tableS1_datasets.csv": table_s1(lines),
        "tableS2_ranking154.csv": table_s2(lines),
        "tableS4_modules.csv": table_s4(lines),
        "tableS6_complexity.csv": table_s6(lines),
        "tableS7_core_signature.csv": table_s7(lines),
        "tableS8_origins.csv": table_s8(lines),
    }
    for name, df in tables.items():
        df.to_csv(OUT / name, index=False)
        print(f"  {name:30s} {len(df):4d} rows  (checks passed)")
    s9b = table_s9b(lines)
    (CATALOG_DATA / "module_criteria_s9b.json").write_text(
        json.dumps({"_source": "Spiro (2026), Supporting Information, Table S9b, as printed (gene lists and "
                               "free-text cells are abbreviated in the print; the complete table is "
                               "tableS9b_module_criteria.csv in the Zenodo record). 'recorded at catalog "
                               "construction' modules have a blank receptor field in the paper.",
                    "modules": dict(sorted(s9b.items()))}, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"  {'module_criteria_s9b.json':30s} {len(s9b):4d} rows  (checks passed)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
