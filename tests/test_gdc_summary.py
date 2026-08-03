"""GDC summary lookup: row->Finding shape + tiering (mini fixture, no download)."""
import sqlite3
from pathlib import Path
from methylask.providers.gdc import GdcProvider


def _mini_summary(tmp_path) -> Path:
    db = tmp_path / "gdc.db"
    con = sqlite3.connect(str(db))
    con.execute("""CREATE TABLE cpg_summary(cpg TEXT, project TEXT, n_tumor INT,
        mean_tumor REAL, n_normal INT, mean_normal REAL)""")
    con.executemany("INSERT INTO cpg_summary VALUES (?,?,?,?,?,?)", [
        # big, well-sampled delta -> robust
        ("cg111", "TCGA-BRCA", 50, 0.80, 10, 0.20),
        # moderate delta, decent tumor n -> moderate
        ("cg222", "TCGA-LUAD", 15, 0.55, 6, 0.42),
        # tiny delta -> speculative
        ("cg333", "TCGA-COAD", 30, 0.51, 8, 0.49),
        # one arm missing -> skipped entirely
        ("cg444", "TCGA-BRCA", 12, 0.60, 0, None),
    ])
    con.execute("CREATE INDEX idx_cpg ON cpg_summary(cpg)")
    con.commit(); con.close()
    return db


def test_robust_hypermethylated(tmp_path):
    p = GdcProvider(summary_path=str(_mini_summary(tmp_path)))
    fs = p.get("cg111")
    assert len(fs) == 1
    assert fs[0].tier.value == "robust"
    assert "hypermethylated" in fs[0].description
    assert fs[0].detail["topic"] == "cancer" and fs[0].detail["project"] == "TCGA-BRCA"


def test_moderate_and_speculative(tmp_path):
    p = GdcProvider(summary_path=str(_mini_summary(tmp_path)))
    assert p.get("cg222")[0].tier.value == "moderate"
    assert p.get("cg333")[0].tier.value == "speculative"


def test_missing_arm_skipped(tmp_path):
    p = GdcProvider(summary_path=str(_mini_summary(tmp_path)))
    assert p.get("cg444") == []      # no normal arm -> no finding


def test_absent_marker_empty(tmp_path):
    p = GdcProvider(summary_path=str(_mini_summary(tmp_path)))
    assert p.get("cg999") == []


def test_summary_path_resolved_from_env(tmp_path, monkeypatch):
    """A provider built with NO argument must still find the summary via
    GDC_SUMMARY_DB. Every production call site does a bare `GdcProvider()`
    (orchestrate.py, cli.py), so when only refresh() read that variable the read
    path returned nothing no matter how complete the summary on disk was — and
    every test above passing summary_path= explicitly is what hid it."""
    monkeypatch.setenv("GDC_SUMMARY_DB", str(_mini_summary(tmp_path)))
    assert GdcProvider().get("cg111")[0].tier.value == "robust"


def test_no_summary_configured_is_empty(monkeypatch):
    """With neither an argument nor the env var, get() stays empty rather than
    raising — status() is what explains the missing summary."""
    monkeypatch.delenv("GDC_SUMMARY_DB", raising=False)
    assert GdcProvider().get("cg111") == []
