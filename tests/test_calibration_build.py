"""Reading the HPA table and building the single_cell_scaled calibration (script 13)."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import perceptome as pct

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "13_build_single_cell_calibration.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("build_sc_calibration", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_hpa(path, cell_types, n_filler=300, seed=0, header=True):
    rng = np.random.default_rng(seed)
    genes = sorted({g for m in pct.load_catalog()["modules"].values() for g in m["core_genes"]})
    genes += [f"FILL{i}" for i in range(n_filler)]
    lines = ["Gene\tGene name\tCell type\tnCPM"] if header else []
    for i, g in enumerate(genes):
        for ct in cell_types:
            lines.append(f"ENSG{i:08d}\t{g}\t{ct}\t{rng.gamma(0.8, 30):.4f}")
    Path(path).write_text("\n".join(lines) + "\n")


def test_reader_returns_expression_and_total_over_all_genes(tmp_path):
    f = tmp_path / "hpa.tsv"
    _write_hpa(f, ["ct_a", "ct_b"])
    expr, lib = pct.read_hpa_single_cell_type(f, genes=["IRF3", "CGAS"])
    assert list(expr.index) == ["IRF3"] and list(expr.columns) == ["ct_a", "ct_b"]
    all_expr, lib_all = pct.read_hpa_single_cell_type(f)
    np.testing.assert_allclose(lib.values, lib_all.values)
    np.testing.assert_allclose(lib_all.values, all_expr.sum(axis=0).values)
    assert (lib > expr.sum(axis=0)).all()          # the library size counts genes the panel does not


def test_reader_handles_headerless_names_and_bad_columns(tmp_path):
    f = tmp_path / "noheader.tsv"
    f.write_text("a\tb\tc\td\nENSG1\tIRF3\tct\t5.0\nENSG2\tIFNB1\tct\t7.0\n")
    expr, lib = pct.read_hpa_single_cell_type(f)             # four columns: positional fallback
    assert expr.loc["IRF3", "ct"] == 5.0 and lib["ct"] == 12.0
    g = tmp_path / "bad.tsv"
    g.write_text("x\ty\n1\t2\n")
    with pytest.raises(ValueError, match="Expected columns"):
        pct.read_hpa_single_cell_type(g)


def test_reader_last_row_wins_for_duplicate_symbols(tmp_path):
    f = tmp_path / "dup.tsv"
    f.write_text("Gene\tGene name\tCell type\tnCPM\nE1\tIRF3\tct\t1.0\nE2\tIRF3\tct\t4.0\n")
    expr, lib = pct.read_hpa_single_cell_type(f)
    assert expr.loc["IRF3", "ct"] == 4.0 and lib["ct"] == 5.0


def test_script_refuses_to_publish_a_calibration_that_does_not_match_table_s7(tmp_path, capsys):
    ref_types = list(pct.load_hpa_perceptivity()["R"].index)
    hpa = tmp_path / "hpa.tsv"
    _write_hpa(hpa, ref_types, n_filler=100)
    out = tmp_path / "calibration_single_cell_scaled_v03.json"
    code = _load_script().main([str(hpa), "--out", str(out)])
    assert code == 2
    assert not out.exists()
    unverified = tmp_path / "calibration_single_cell_scaled_v03.UNVERIFIED.json"
    assert unverified.exists()
    cal = pct.load_calibration(str(unverified))
    assert cal.name == "single_cell_scaled" and cal.reference_scores.shape == (154, 44)
    assert cal.centroids.shape == (154, 9)
    assert "Does NOT reproduce Table S7" in capsys.readouterr().out


def test_script_publishes_when_the_table_s7_check_passes(tmp_path, monkeypatch):
    ref_types = list(pct.load_hpa_perceptivity()["R"].index)
    hpa = tmp_path / "hpa.tsv"
    _write_hpa(hpa, ref_types, n_filler=100)
    script = _load_script()
    monkeypatch.setattr(script, "check_calibration_against_table_s7",
                        lambda cal: {"megakaryocyte_r": 1.0, "cytotrophoblast_r": 1.0,
                                     "max_abs_diff": 0.0, "mean_abs_diff": 0.0, "max_abs_diff_checked": 0.0,
                                     "n_checked": 86, "n_exact": 86, "known_deviations": {}, "matches": True})
    out = tmp_path / "cal.json"
    assert script.main([str(hpa), "--out", str(out)]) == 0
    assert out.exists()
    assert pct.load_calibration(str(out)).name == "single_cell_scaled"


def test_script_stops_when_cell_types_are_missing(tmp_path):
    hpa = tmp_path / "hpa.tsv"
    _write_hpa(hpa, ["only_one_cell_type"], n_filler=10)
    with pytest.raises(SystemExit, match="154 reference cell types"):
        _load_script().main([str(hpa), "--out", str(tmp_path / "x.json")])
