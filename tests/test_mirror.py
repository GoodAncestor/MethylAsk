"""EWAS mirror: mirror-first lookup + row->Finding shape (mini fixture DB)."""
import sqlite3, os
from pathlib import Path
from methylask.providers import ewas_mirror as M


def _mini_db(tmp_path) -> Path:
    db = tmp_path / "mini.db"
    con = sqlite3.connect(str(db))
    con.execute("""CREATE TABLE findings(cpg TEXT, trait TEXT, gene TEXT, beta TEXT,
        se TEXT, p TEXT, n TEXT, tissue TEXT, methylation_array TEXT, chrpos TEXT,
        pmid TEXT, efo TEXT)""")
    con.execute("INSERT INTO findings VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                ("cg99999999", "age", "TESTG", "0.1", "0.01", "1e-9", "5000",
                 "Whole blood", "EPIC", "chr1:1000", "12345678", "EFO_0000246"))
    con.execute("CREATE INDEX idx_cpg ON findings(cpg)")
    con.commit(); con.close()
    return db


def test_mirror_lookup_hit(tmp_path):
    db = _mini_db(tmp_path)
    rows = M.mirror_lookup("cg99999999", db)
    assert len(rows) == 1
    assert rows[0]["trait"] == "age" and rows[0]["n"] == "5000"
    assert rows[0]["efo"] == "EFO_0000246"


def test_mirror_lookup_miss_returns_empty_not_none(tmp_path):
    db = _mini_db(tmp_path)
    # a CpG absent from the mirror returns [] (mirror exists), not None
    assert M.mirror_lookup("cg00000000", db) == []


def test_mirror_absent_returns_none(tmp_path):
    # no db at all -> None, so the provider falls back to cache/live
    assert M.mirror_lookup("cg99999999", tmp_path / "nope.db") is None


def test_provider_uses_mirror(tmp_path, monkeypatch):
    db = _mini_db(tmp_path)
    monkeypatch.setattr(M, "MIRROR_DB", db)
    from methylask.providers.ewas_catalog import EwasCatalogProvider
    fs = EwasCatalogProvider().get("cg99999999")
    assert len(fs) == 1
    assert fs[0].tier.value == "robust"          # n=5000, p=1e-9
    assert fs[0].detail["topic"] == "aging"
    assert fs[0].detail["efo"] == "EFO_0000246"
