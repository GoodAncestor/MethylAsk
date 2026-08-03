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


def test_stratified_plan_balances_arms_and_drops_one_armed_projects(monkeypatch):
    """The sampler exists because GDC's page order is clustered by submission
    batch, not shuffled — the first 500 files of the real corpus are two projects.
    So it must (a) reach every project, (b) cap each ARM separately rather than
    each project, and (c) drop projects with no normal arm, whose files can only
    cost bandwidth: _get_from_summary skips any cell missing an arm."""
    corpus = (
        [{"file_id": f"t-brca-{i}", "project": "TCGA-BRCA", "is_normal": False} for i in range(100)]
        + [{"file_id": f"n-brca-{i}", "project": "TCGA-BRCA", "is_normal": True} for i in range(5)]
        # last in page order, and thinly represented — a first-N slice would miss it
        + [{"file_id": "t-chol-0", "project": "TCGA-CHOL", "is_normal": False},
           {"file_id": "n-chol-0", "project": "TCGA-CHOL", "is_normal": True}]
        # tumour only: 18 of the 46 real projects look like this
        + [{"file_id": f"t-lgg-{i}", "project": "TCGA-LGG", "is_normal": False} for i in range(50)]
    )
    p = GdcProvider(summary_path="/nonexistent")
    monkeypatch.setattr(p, "_list_files", lambda max_files: iter(corpus))
    chosen = p._plan_stratified(per_arm=3)
    by_project = {}
    for f in chosen:
        by_project.setdefault(f["project"], []).append(f)
    assert set(by_project) == {"TCGA-BRCA", "TCGA-CHOL"}      # LGG dropped, CHOL kept
    assert sum(1 for f in by_project["TCGA-BRCA"] if not f["is_normal"]) == 3
    assert sum(1 for f in by_project["TCGA-BRCA"] if f["is_normal"]) == 3
    assert len(by_project["TCGA-CHOL"]) == 2                  # cap is a ceiling, not a quota


def test_findings_carry_the_sampling_provenance(tmp_path):
    """A sampled summary and a full one are the same file on disk. Without this the
    per-arm counts are the only clue, and they read as 'this CpG is thinly
    measured' rather than 'the whole summary is a sample'."""
    db = _mini_summary(tmp_path)
    con = sqlite3.connect(str(db))
    con.execute("CREATE TABLE summary_meta(key TEXT PRIMARY KEY, value TEXT)")
    con.execute("INSERT INTO summary_meta VALUES ('sampling','stratified: up to 50 files per arm per project')")
    con.commit(); con.close()
    f = GdcProvider(summary_path=str(db)).get("cg111")[0]
    assert "stratified" in f.detail["sampling"]


def test_summary_without_provenance_still_reads(tmp_path):
    """Summaries built before summary_meta existed must keep working — the key is
    simply absent from detail rather than the read blowing up."""
    f = GdcProvider(summary_path=str(_mini_summary(tmp_path))).get("cg111")[0]
    assert "sampling" not in f.detail


def test_no_summary_configured_is_empty(monkeypatch):
    """With neither an argument nor the env var, get() stays empty rather than
    raising — status() is what explains the missing summary."""
    monkeypatch.delenv("GDC_SUMMARY_DB", raising=False)
    assert GdcProvider().get("cg111") == []


def test_unwritable_summary_path_fails_before_downloading(tmp_path, monkeypatch):
    """The writability check must run BEFORE any download. A build dir created by
    docker as root under a container running as the service account used to
    surface as `unable to open database file` only at the very end — after the
    whole corpus had been streamed."""
    called = []
    p = GdcProvider(summary_path=str(tmp_path / "ro" / "gdc.db"))
    (tmp_path / "ro").mkdir(mode=0o500)
    monkeypatch.setattr(p, "_list_files", lambda *a, **k: called.append(1) or iter([]))
    st = p.refresh(per_arm=1)
    assert st.health.value == "unavailable" and "not writable" in st.note
    assert not called       # nothing was listed, so nothing was downloaded
