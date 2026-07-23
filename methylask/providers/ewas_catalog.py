"""EWAS Catalog provider (MRC-IEU).

Validated live 2026-07-23: GET https://www.ewascatalog.org/api/?cpg=<probe>
returns {"fields": [...], "results": [[...], ...]} where each result row is a
POSITIONAL array aligned to `fields` (not a dict) — see docs/VALIDATION.md.

Design posture (docs/DESIGN.md §3.4): the full results dump (174 MB) is mirrored
to local disk by refresh(); per-CpG lookups then hit the local copy. The live
API path here is used for prototyping and as a fallback before the mirror exists.
"""
from __future__ import annotations
import json, ssl, urllib.request, urllib.error
from biocore.providers.base import Provider, Finding, Tier, Category, ProviderStatus, Health

_API = "https://www.ewascatalog.org/api/?cpg="
# Prototype only: some MRC-IEU hosts have intermittent cert issues. Public
# read-only reference data, no user values on the wire. Off by default.
_INSECURE = False


def _tier_from(row: dict) -> Tier:
    """Map study metadata to an evidence tier (docs/DESIGN.md §4.3.1)."""
    try:
        n = int(float(row.get("n") or 0))
    except (TypeError, ValueError):
        n = 0
    try:
        p = float(row.get("p"))
    except (TypeError, ValueError):
        p = 1.0
    if n >= 1000 and p <= 1e-7:
        return Tier.ROBUST
    if n >= 200 and p <= 1e-4:
        return Tier.MODERATE
    return Tier.SPECULATIVE


class EwasCatalogProvider(Provider):
    name = "ewas_catalog"

    def __init__(self, insecure: bool = _INSECURE, timeout: int = 30):
        self._ctx = ssl._create_unverified_context() if insecure else None
        self._timeout = timeout

    def _fetch(self, cpg: str) -> dict:
        req = urllib.request.Request(_API + cpg,
            headers={"User-Agent": "methylask", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=self._timeout, context=self._ctx) as r:
            return json.loads(r.read().decode("utf-8", "replace"))

    def get(self, marker: str) -> list[Finding]:
        try:
            payload = self._fetch(marker)
        except Exception:
            return []  # errors surface via status(), not as exceptions here
        fields = payload.get("fields") or []
        rows = payload.get("results") or []
        out: list[Finding] = []
        for raw in rows:
            row = dict(zip(fields, raw))
            trait = row.get("trait", "unknown trait")
            gene = row.get("gene") or "?"
            out.append(Finding(
                marker=marker, source=self.name,
                description=f"Associated with '{trait}' (gene {gene})",
                tier=_tier_from(row),
                categories=[Category.CLINICAL, Category.TRAIT],
                detail={k: row.get(k) for k in
                        ("beta", "se", "p", "n", "tissue", "methylation_array", "chrpos")},
                link=f"https://www.ewascatalog.org/?query={marker}",
                pmids=[str(row["pmid"])] if row.get("pmid") else [],
            ))
        return out

    def refresh(self) -> ProviderStatus:
        # TODO: download ewascatalog-results.txt.gz (174 MB) -> local SQLite/Parquet
        return self.status()

    def status(self) -> ProviderStatus:
        try:
            self._fetch("cg00000029")
            return ProviderStatus(self.name, Health.OK, note="live API reachable")
        except urllib.error.HTTPError as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE, note=f"HTTP {e.code}")
        except Exception as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE,
                                  note=f"{type(e).__name__}: {str(e)[:80]}")
