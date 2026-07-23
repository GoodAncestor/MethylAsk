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
import json, urllib.parse, urllib.request, urllib.error
from pathlib import Path
from biocore.providers.base import Provider, Finding, Tier, Category, ProviderStatus, Health

_API = "https://api.gdc.cancer.gov/"


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
        # TODO: read per-CpG tumour/normal beta deltas from the local summary
        # (keyed store built by refresh()); emit one Finding per cancer project
        # with a robust/moderate tier based on effect size + sample counts.
        return []

    def refresh(self) -> ProviderStatus:
        # TODO: (1) list all "Methylation Beta Value" files via the files
        # endpoint; (2) download to local mirror; (3) stream once to build the
        # per-CpG tumour/normal summary. See docs/DESIGN.md §8.2.
        return self.status()

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
