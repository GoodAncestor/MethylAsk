"""EWAS Catalog provider (MRC-IEU).

Validated live 2026-07-23: GET https://www.ewascatalog.org/api/?cpg=<probe>
returns {"fields": [...], "results": [[...], ...]} where each result row is a
POSITIONAL array aligned to `fields` (not a dict) — see docs/VALIDATION.md.

Design posture (docs/DESIGN.md §3.4): the full results dump (174 MB) is mirrored
to local disk by refresh(); per-CpG lookups then hit the local copy. The live
API path here is used for prototyping and as a fallback before the mirror exists.
"""
from __future__ import annotations
import json, ssl, os, hashlib, time, urllib.request, urllib.error
from pathlib import Path
from biocore.providers.base import Provider, Finding, Tier, Category, ProviderStatus, Health
from ..traits import classify_topic, humanize_trait

# per-CpG response cache on disk. Live EWAS lookups are ~1 HTTP call per marker
# (slow: 10-13s for a 40-marker report). The API response for a given CpG is
# effectively static between catalog releases, so cache it. Dir is configurable
# (DNAREPORT/METHYLASK_CACHE_DIR) and defaults to a temp path; TTL is generous.
_CACHE_DIR = Path(os.environ.get("METHYLASK_CACHE_DIR")
                  or os.environ.get("DNAREPORT_CACHE_DIR")
                  or (Path(os.environ.get("TMPDIR", "/tmp")) / "methylask_ewas_cache"))
_CACHE_TTL = int(os.environ.get("METHYLASK_CACHE_TTL", str(30 * 24 * 3600)))  # 30 days

# map a fine-grained topic to the coarse Category used for report sections.
# aging -> AGING; disease/biomarker topics -> CLINICAL; the rest -> TRAIT.
_TOPIC_CATEGORY = {
    "aging": Category.AGING,
    "cancer": Category.CLINICAL,
    "cardiovascular": Category.CLINICAL,
    "immune": Category.CLINICAL,
    "respiratory": Category.CLINICAL,
    "neuro": Category.CLINICAL,
    "metabolic": Category.CLINICAL,
    "reproductive": Category.TRAIT,
    "lifestyle": Category.TRAIT,
    "proteomic": Category.TRAIT,
    "other": Category.TRAIT,
}

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

    def _cache_file(self, cpg: str) -> Path:
        h = hashlib.sha1(cpg.encode()).hexdigest()[:16]
        return _CACHE_DIR / f"{h}.json"

    def _fetch(self, cpg: str) -> dict:
        # disk cache: return a fresh cached response if present
        cf = self._cache_file(cpg)
        try:
            if cf.exists() and (time.time() - cf.stat().st_mtime) < _CACHE_TTL:
                return json.loads(cf.read_text())
        except Exception:
            pass  # corrupt/unreadable cache -> fall through to live fetch
        req = urllib.request.Request(_API + cpg,
            headers={"User-Agent": "methylask", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=self._timeout, context=self._ctx) as r:
            payload = json.loads(r.read().decode("utf-8", "replace"))
        # write-through (best-effort; never fail a lookup on a cache write error)
        try:
            _CACHE_DIR.mkdir(parents=True, exist_ok=True)
            cf.write_text(json.dumps(payload))
        except Exception:
            pass
        return payload

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
            gene = (row.get("gene") or "").strip() or None

            topic = classify_topic(trait)
            label, kind, accession = humanize_trait(trait)
            # a plain-language description; a protein-accession trait reads as a
            # "protein level" measurement rather than a bare code
            if kind == "protein":
                desc = f"linked to blood level of protein {label}"
            else:
                desc = f"linked to {label}"

            # coarse Category for section placement (renderer groups by these);
            # the fine-grained `topic` drives the subject filter.
            cats = [_TOPIC_CATEGORY.get(topic, Category.TRAIT)]

            detail = {k: row.get(k) for k in
                      ("beta", "se", "p", "n", "tissue", "methylation_array", "chrpos")}
            detail["topic"] = topic
            detail["trait"] = label
            if gene:
                detail["gene"] = gene
            if accession:
                detail["protein"] = accession

            out.append(Finding(
                marker=marker, source=self.name,
                description=desc,
                tier=_tier_from(row),
                categories=cats,
                detail=detail,
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
