"""GDC / TCGA cancer-methylation provider (NCI Genomic Data Commons).

Validated live 2026-07-23: api.gdc.cancer.gov returns 20,397 "Methylation Beta
Value" files (SeSAMe level-3 beta, EPIC/450K) across ~33 cancer projects.

Two-layer design (docs/DESIGN.md §8.2):
  - refresh() streams the corpus once (20,397 files, 293.8 GB measured 2026-08-03)
    and builds a per-CpG summary: for each probe, the tumour-vs-normal beta
    distribution per cancer project. Output is a few hundred MB and serves fast
    lookups forever.
  - get() reads that local summary. Before the summary exists, get() returns
    nothing and status() says why.

There is no per-CpG endpoint upstream — the only calls are `files` (metadata),
`data/<file_id>` (a whole per-sample beta matrix) and `status`. Answering one CpG
means reading the corpus, which is why the summary is built ahead of time rather
than fetched on demand.

**Sampling.** `per_arm=N` builds from a stratified slice instead of all 293.8 GB:
up to N tumour and N normal files from every project that HAS both arms. That is
hours instead of days and covers the same 28 projects, because the corpus is far
thinner than its file count suggests — 18 of the 46 projects have NO normal arm at
all, so they can never yield a tumour-vs-normal answer no matter how much is
downloaded, and the normal arms that do exist run from 612 files (CPTAC-3) down to
2 (TCGA-SKCM, GBM, THYM). Measured 2026-08-03: N=50 is 2,185 files / 25.1 GB.

Do NOT use `max_files=N` for this. It takes the first N files in API page order,
and that order is clustered by submission batch, not shuffled — the first 500 files
are two projects (CPTAC-3, HCMI-CMDC), and HCMI-CMDC has no normals. `max_files`
exists for validation runs; it is a slice, not a sample.

A sampled build records what it did in `summary_meta` and every Finding it produces
carries a `sampling` key, so a thin summary can't be mistaken for a full one.
"""
from __future__ import annotations
import os, json, math, sqlite3, tempfile, time, urllib.parse, urllib.request, urllib.error
from pathlib import Path
from biocore.providers.base import Provider, Finding, Tier, Category, ProviderStatus, Health

_API = "https://api.gdc.cancer.gov/"
# sample_type strings GDC uses; anything containing "Normal" is the normal arm.
# Substring, so it also catches "bone marrow normal" / "blood derived normal" /
# "fibroblasts from bone marrow normal". 2,839 of the 20,397 files carry NO
# sample_type at all; those fall to the tumour arm, which is the conservative way
# round — an untyped file inflating the tumour n cannot invent a normal arm.
_NORMAL_HINT = "normal"
# Conventional in-container path for the summary volume, used only when neither a
# constructor argument nor GDC_SUMMARY_DB says otherwise.
_DEFAULT_SUMMARY_PATH = "/data/gdc_summary/gdc_summary.db"


class GdcProvider(Provider):
    name = "gdc"

    # Where the prebuilt summary lives when no path is passed. A container mounts
    # the summary volume and sets this; refresh() has always honoured it, so the
    # READ path must resolve it identically or the two disagree.
    _SUMMARY_ENV = "GDC_SUMMARY_DB"

    def __init__(self, summary_path: str | None = None, timeout: int = 60):
        # Fall back to the environment. Every production call site constructs a
        # bare GdcProvider() (dnareport/orchestrate.py, methylask/cli.py), so
        # while only refresh() read GDC_SUMMARY_DB, get() returned [] no matter
        # how complete the summary on disk was — a built mirror was unreachable.
        resolved = summary_path or os.environ.get(self._SUMMARY_ENV)
        self._summary_path = Path(resolved) if resolved else None
        self._timeout = timeout
        self._meta_cache: dict | None = None

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
        sampling = self._meta().get("sampling")
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
                        "mean_normal": round(r["mean_normal"], 4),
                        # what the summary underneath this number was built from,
                        # so a reader can weigh it without going and looking
                        **({"sampling": sampling} if sampling else {})},
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

    def _plan_stratified(self, per_arm: int) -> list[dict]:
        """Pick up to `per_arm` tumour and `per_arm` normal files from every project
        that has both arms, and return them as a concrete list.

        Listing the whole corpus first costs ~41 metadata pages and no beta data,
        which is what buys an even sample: taking the first N files instead would
        follow GDC's page order, and that order is clustered by submission batch.

        Projects with an empty normal arm are dropped here rather than downloaded
        and discarded later — _get_from_summary() skips any cell missing an arm, so
        those files can only ever cost bandwidth."""
        by: dict[str, dict[str, list]] = {}
        for f in self._list_files(None):
            arms = by.setdefault(f["project"], {"t": [], "n": []})
            arms["n" if f["is_normal"] else "t"].append(f)
        chosen: list[dict] = []
        for proj in sorted(by):
            arms = by[proj]
            if not arms["t"] or not arms["n"]:
                continue
            chosen.extend(arms["t"][:per_arm])
            chosen.extend(arms["n"][:per_arm])
        return chosen

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

    def refresh(self, max_files: int | None = None, workdir: str | None = None,
                per_arm: int | None = None) -> ProviderStatus:
        """Build the per-CpG tumour/normal summary. Streams each beta file once,
        accumulating per (cpg, project, arm) running count+mean, into a SQLite at
        self._summary_path. Runs on a worker via `worker.py refresh:gdc`.

        Three ways to choose the input, in order of preference:
          per_arm=N    stratified sample, N per arm per project — the useful one
          max_files=N  the first N files in page order — validation only, biased
          neither      the whole 293.8 GB corpus

        Peak memory is the accumulator, one entry per (cpg, project) cell: ~865k
        EPIC probes × 28 projects for a stratified build, and it is the reason
        refresh-publish.sh caps the container at 24 GB."""
        # __init__ already resolved an argument or GDC_SUMMARY_DB; reaching here
        # means neither was given. Don't re-read the env — that is how the read and
        # write paths drifted apart in the first place.
        if not self._summary_path:
            self._summary_path = Path(_DEFAULT_SUMMARY_PATH)
        wd = Path(workdir or tempfile.mkdtemp(prefix="gdc-"))
        wd.mkdir(parents=True, exist_ok=True)
        if per_arm:
            plan = self._plan_stratified(per_arm)
            sampling = f"stratified: up to {per_arm} files per arm per project"
        elif max_files:
            plan = self._list_files(max_files)
            sampling = f"first {max_files} files in API page order (biased — validation only)"
        else:
            plan = self._list_files(None)
            sampling = "full corpus"
        # accumulator: {(cpg, project): [n_t, sum_t, n_n, sum_n]}
        acc: dict[tuple, list] = {}
        n_files = 0
        try:
            for f in plan:
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
                local.unlink(missing_ok=True)   # stream: don't keep the full corpus
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
        # Provenance. Without it a 2,185-file sample and a 20,397-file full build
        # are the same file on disk, and the per-arm counts in each row are the
        # only clue — which reads as "this CpG happens to be thinly measured"
        # rather than "the whole summary is a sample". status() and every Finding
        # read this back.
        projects = sorted({proj for _, proj in acc})
        con.execute("DROP TABLE IF EXISTS summary_meta")
        con.execute("CREATE TABLE summary_meta(key TEXT PRIMARY KEY, value TEXT)")
        con.executemany("INSERT INTO summary_meta VALUES (?,?)", [
            ("sampling", sampling),
            ("n_files", str(n_files)),
            ("n_projects", str(len(projects))),
            ("projects", ",".join(projects)),
            ("built_at", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
            ("data_release", self._data_release() or "?"),
        ])
        con.commit(); con.close()
        return ProviderStatus(self.name, Health.OK,
            note=f"summary built: {len(acc)} cpg×project cells from {n_files} files "
                 f"across {len(projects)} projects ({sampling})")

    def _data_release(self) -> str | None:
        try:
            return self._api("status", {}).get("data_release")
        except Exception:
            return None   # provenance is best-effort; never fail a finished build

    def _meta(self) -> dict:
        """summary_meta as a dict, read once and cached. Absent table -> {}, which
        is what a summary built before provenance existed looks like."""
        if self.__dict__.get("_meta_cache") is None:
            m = {}
            if self._summary_ready():
                con = sqlite3.connect(str(self._summary_path))
                try:
                    m = dict(con.execute("SELECT key, value FROM summary_meta"))
                except sqlite3.OperationalError:
                    m = {}
                finally:
                    con.close()
            self._meta_cache = m
        return self._meta_cache

    def status(self) -> ProviderStatus:
        try:
            st = self._api("status", {})
            rel = st.get("data_release", "?")
            if self._summary_ready():
                m = self._meta()
                built = f"{m.get('n_files', '?')} files / {m.get('n_projects', '?')} projects"
                what = m.get("sampling", "provenance not recorded (pre-2026-08 build)")
                return ProviderStatus(
                    self.name, Health.OK, version=rel,
                    record_count=int(m["n_files"]) if m.get("n_files", "").isdigit() else None,
                    note=f"local CpG summary ready — {built}, {what}")
            return ProviderStatus(self.name, Health.STALE, version=rel,
                                  note="API reachable; per-CpG summary not built (run refresh)")
        except urllib.error.HTTPError as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE, note=f"HTTP {e.code}")
        except Exception as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE,
                                  note=f"{type(e).__name__}: {str(e)[:80]}")
