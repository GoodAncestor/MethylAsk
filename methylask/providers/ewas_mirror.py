"""Local EWAS Catalog mirror — build once, query offline.

The live API is one HTTP call per CpG (slow, and a dependency on MRC-IEU uptime).
This mirrors the two bulk downloads into a local SQLite keyed by CpG so lookups
are instant and offline:

  - ewascatalog-results.txt.gz : per-association rows
      cpg, beta, se, p, details, study_id, cpg, chrpos, chr, pos, gene, type, assocs
  - ewascatalog-studies.txt.gz : per-study metadata (joined on study_id)
      author, pmid, trait, efo, ..., methylation_array, tissue, ..., n, ...

Joined on study_id, this reconstructs the same fields the API returns
(trait, gene, beta, se, p, n, tissue, methylation_array, chrpos, pmid) plus the
EFO ontology id. build_mirror() streams both files so peak memory stays bounded.
"""
from __future__ import annotations
import os, gzip, sqlite3, urllib.request
from pathlib import Path

_BASE = "https://www.ewascatalog.org/static/docs"
_RESULTS_URL = f"{_BASE}/ewascatalog-results.txt.gz"
_STUDIES_URL = f"{_BASE}/ewascatalog-studies.txt.gz"

# where the mirror db lives; configurable so a worker/NAS can host a shared copy
MIRROR_DB = Path(os.environ.get("EWAS_MIRROR_DB")
                 or (Path(os.environ.get("METHYLASK_DATA", "/tmp")) / "ewas_mirror.db"))


def _download(url: str, dest: Path):
    req = urllib.request.Request(url, headers={"User-Agent": "methylask"})
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as fh:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            fh.write(b)


def _rows(path: Path):
    """Yield dict rows from a gzipped TSV (header-driven)."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            vals = line.rstrip("\n").split("\t")
            if len(vals) >= len(header):
                yield dict(zip(header, vals))


def build_mirror(db_path: Path | None = None, workdir: Path | None = None) -> dict:
    """Download both bulk files and build the SQLite mirror. Returns a summary
    {n_studies, n_findings, db_path}. Idempotent: rebuilds the table each run."""
    db_path = Path(db_path or MIRROR_DB)
    workdir = Path(workdir or db_path.parent)
    workdir.mkdir(parents=True, exist_ok=True)
    results_gz = workdir / "ewascatalog-results.txt.gz"
    studies_gz = workdir / "ewascatalog-studies.txt.gz"
    _download(_RESULTS_URL, results_gz)
    _download(_STUDIES_URL, studies_gz)

    # 1. load study metadata into a dict keyed by study_id (small: ~thousands)
    studies: dict[str, dict] = {}
    for s in _rows(studies_gz):
        studies[s.get("study_id", "")] = s

    # 2. stream results, join study metadata, write to SQLite keyed by cpg
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(db_path))
    con.execute("DROP TABLE IF EXISTS findings")
    con.execute("""CREATE TABLE findings(
        cpg TEXT, trait TEXT, gene TEXT, beta TEXT, se TEXT, p TEXT, n TEXT,
        tissue TEXT, methylation_array TEXT, chrpos TEXT, pmid TEXT, efo TEXT)""")
    n = 0
    batch = []
    for r in _rows(results_gz):
        st = studies.get(r.get("study_id", ""), {})
        batch.append((
            r.get("cpg"), st.get("trait", "unknown trait"), r.get("gene"),
            r.get("beta"), r.get("se"), r.get("p"), st.get("n"),
            st.get("tissue"), st.get("methylation_array"), r.get("chrpos"),
            st.get("pmid"), st.get("efo")))
        if len(batch) >= 5000:
            con.executemany("INSERT INTO findings VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", batch)
            n += len(batch); batch = []
    if batch:
        con.executemany("INSERT INTO findings VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", batch)
        n += len(batch)
    con.execute("CREATE INDEX idx_cpg ON findings(cpg)")
    con.commit(); con.close()
    return {"n_studies": len(studies), "n_findings": n, "db_path": str(db_path)}


def mirror_lookup(cpg: str, db_path: Path | None = None) -> list[dict] | None:
    """Return per-association dict rows for a CpG from the mirror, or None if no
    mirror exists (caller then falls back to cache/live)."""
    db_path = Path(db_path or MIRROR_DB)
    if not db_path.exists():
        return None
    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(
            "SELECT trait,gene,beta,se,p,n,tissue,methylation_array,chrpos,pmid,efo "
            "FROM findings WHERE cpg=?", (cpg,)).fetchall()
    except sqlite3.OperationalError:
        return None
    finally:
        con.close()
    return [dict(r) for r in rows]
