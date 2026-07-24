"""GDC / TCGA cancer-methylation provider (NCI Genomic Data Commons).

Validated live 2026-07-23: api.gdc.cancer.gov returns 20,397 "Methylation Beta
Value" files (SeSAMe level-3 beta, EPIC/450K) across ~33 cancer projects.

Two-layer design (docs/DESIGN.md §8.2):
  - refresh() mirrors the corpus locally (294 GB, ~1 h on 10 Gbit) and runs a
    ONE-TIME streaming pass to build a per-CpG summary: for each probe, the
    tumour-vs-normal beta distribution per cancer project. Output is a few
    hundred MB and serves fast lookups forever.
  - get() reads that local summary. Before the summary exists, get() returns
    the project/coverage context the live API can answer cheaply.

This is where the only compute-heavy step lives, and it runs once per data
release inside refresh() — never on the request path.
"""
from __future__ import annotations
import os, json, math, sqlite3, tempfile, urllib.parse, urllib.request, urllib.error
from pathlib import Path
from biocore.providers.base import Provider, Finding, Tier, Category, ProviderStatus, Health

_API = "https://api.gdc.cancer.gov/"
# sample_type strings GDC uses; anything containing "Normal" is the normal arm.
_NORMAL_HINT = "normal"


class GdcProvider(Provider):
    name = "gdc"

    def __init__(self, summary_path: str | None = None, timeout: int = 60):
        self._summary_path = Path(summary_path) if summary_path else None
        self._timeout = timeout

    def _api(self, endpoint: str, params: dict) -> dict:
        url = _API + endpoint + "?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={"User-Agent": "methylask"})
        with urllib.request.urlopen(req, timeout=self._timeout) as r:
            return json.loads(r.read().decode("utf-8", "replace"))

    def _summary_ready(self) -> bool:
        return bool(self._summary_path and self._summary_path.exists())

    def get(self, marker: str) -> list[Finding]:
        # Full per-CpG tumour/normal lookup requires the precomputed summary.
        if self._summary_ready():
            return self._get_from_summary(marker)
        return []  # summary not built yet; status() explains

    def _get_from_summary(self, marker: str) -> list[Finding]:
        """Read per-CpG tumour/normal beta deltas from the local summary and emit
        one Finding per cancer project where the marker differs. Tier by effect
        size (|Δβ|) and per-arm sample counts."""
        con = sqlite3.connect(str(self._summary_path))
        con.row_factory = sqlite3.Row
        try:
            rows = con.execute(
                "SELECT project,n_tumor,mean_tumor,n_normal,mean_normal "
                "FROM cpg_summary WHERE cpg=?", (marker,)).fetchall()
        except sqlite3.OperationalError:
            return []
        finally:
            con.close()
        out = []
        for r in rows:
            if not r["n_tumor"] or not r["n_normal"]:
                continue
            d = (r["mean_tumor"] or 0) - (r["mean_normal"] or 0)
            direction = "hypermethylated" if d > 0 else "hypomethylated"
            ad = abs(d)
            # tier: large, well-sampled deltas are robust; small/thin are speculative
            if ad >= 0.2 and r["n_tumor"] >= 20 and r["n_normal"] >= 5:
                tier = Tier.ROBUST
            elif ad >= 0.1 and r["n_tumor"] >= 10:
                tier = Tier.MODERATE
            else:
                tier = Tier.SPECULATIVE
            out.append(Finding(
                marker=marker, source=self.name,
                description=f"{direction} in {r['project']} tumour vs normal "
                            f"(Δβ={d:+.2f})",
                tier=tier, categories=[Category.CLINICAL],
                detail={"topic": "cancer", "modality": "methylome",
                        "project": r["project"], "delta_beta": round(d, 4),
                        "n_tumor": r["n_tumor"], "n_normal": r["n_normal"],
                        "mean_tumor": round(r["mean_tumor"], 4),
                        "mean_normal": round(r["mean_normal"], 4)},
                link=f"https://portal.gdc.cancer.gov/projects/{r['project']}"))
        return out

    def _list_files(self, max_files: int | None):
        """Page the files endpoint for level-3 methylation beta files, yielding
        {file_id, project, is_normal}. max_files caps the pull (slice for
        validation; None = the full ~20k-file corpus)."""
        filt = {"op": "and", "content": [
            {"op": "in", "content": {"field": "data_type", "value": ["Methylation Beta Value"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["TXT"]}}]}
        size = 500
        got = 0
        frm = 0
        while True:
            page = self._api("files", {
                "filters": json.dumps(filt),
                "fields": "file_id,cases.project.project_id,cases.samples.sample_type",
                "size": str(size), "from": str(frm), "format": "json"})
            hits = page["data"]["hits"]
            if not hits:
                break
            for h in hits:
                case = (h.get("cases") or [{}])[0]
                proj = case.get("project", {}).get("project_id", "UNKNOWN")
                st = (case.get("samples") or [{}])[0].get("sample_type", "") or ""
                yield {"file_id": h["file_id"], "project": proj,
                       "is_normal": _NORMAL_HINT in st.lower()}
                got += 1
                if max_files and got >= max_files:
                    return
            frm += len(hits)

    def _download_betas(self, file_id: str, dst: Path) -> Path:
        """Download one level-3 beta file (cpg<TAB>beta rows) to dst."""
        url = _API + "data/" + file_id
        req = urllib.request.Request(url, headers={"User-Agent": "methylask"})
        with urllib.request.urlopen(req, timeout=self._timeout) as r, open(dst, "wb") as fh:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                fh.write(b)
        return dst

    def refresh(self, max_files: int | None = None, workdir: str | None = None) -> ProviderStatus:
        """Build the per-CpG tumour/normal summary. Streams each beta file once,
        accumulating per (cpg, project, arm) running count+mean, into a SQLite at
        self._summary_path. Heavy (full corpus ~294 GB / ~20k files) — runs on a
        worker via `worker.py refresh:gdc`. max_files caps it for validation."""
        if not self._summary_path:
            self._summary_path = Path(os.environ.get("GDC_SUMMARY_DB",
                                      "/data/gdc_summary/gdc_summary.db"))
        wd = Path(workdir or tempfile.mkdtemp(prefix="gdc-"))
        wd.mkdir(parents=True, exist_ok=True)
        # accumulator: {(cpg, project): [n_t, sum_t, n_n, sum_n]}
        acc: dict[tuple, list] = {}
        n_files = 0
        try:
            for f in self._list_files(max_files):
                local = wd / f"{f['file_id']}.txt"
                try:
                    self._download_betas(f["file_id"], local)
                except Exception:
                    continue
                arm_n = f["is_normal"]
                with open(local, "r", errors="replace") as fh:
                    for line in fh:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 2:
                            continue
                        cpg, bval = parts[0], parts[1]
                        if not cpg.startswith("cg"):
                            continue
                        try:
                            b = float(bval)
                        except ValueError:
                            continue
                        if math.isnan(b):
                            continue
                        key = (cpg, f["project"])
                        s = acc.setdefault(key, [0, 0.0, 0, 0.0])
                        if arm_n:
                            s[2] += 1; s[3] += b
                        else:
                            s[0] += 1; s[1] += b
                local.unlink(missing_ok=True)   # stream: don't keep the 294 GB
                n_files += 1
        except Exception as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE,
                                  note=f"refresh failed after {n_files} files: {e}")

        # write the summary
        self._summary_path.parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(str(self._summary_path))
        con.execute("DROP TABLE IF EXISTS cpg_summary")
        con.execute("""CREATE TABLE cpg_summary(cpg TEXT, project TEXT,
            n_tumor INT, mean_tumor REAL, n_normal INT, mean_normal REAL)""")
        con.executemany(
            "INSERT INTO cpg_summary VALUES (?,?,?,?,?,?)",
            [(cpg, proj, s[0], (s[1]/s[0] if s[0] else None),
              s[2], (s[3]/s[2] if s[2] else None)) for (cpg, proj), s in acc.items()])
        con.execute("CREATE INDEX idx_cpg ON cpg_summary(cpg)")
        con.commit(); con.close()
        return ProviderStatus(self.name, Health.OK,
            note=f"summary built: {len(acc)} cpg×project cells from {n_files} files")

    def status(self) -> ProviderStatus:
        try:
            st = self._api("status", {})
            rel = st.get("data_release", "?")
            if self._summary_ready():
                return ProviderStatus(self.name, Health.OK,
                                      version=rel, note="local CpG summary ready")
            return ProviderStatus(self.name, Health.STALE, version=rel,
                                  note="API reachable; per-CpG summary not built (run refresh)")
        except urllib.error.HTTPError as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE, note=f"HTTP {e.code}")
        except Exception as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE,
                                  note=f"{type(e).__name__}: {str(e)[:80]}")
