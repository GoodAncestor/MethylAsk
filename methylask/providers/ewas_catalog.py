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
from ..traits import classify_topic, humanize_trait, trait_class, trait_copy_key
from .ewas_mirror import mirror_lookup, build_mirror

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

    def _rows_for(self, marker: str) -> list[dict]:
        """Per-association rows for a CpG, mirror-first: a built local mirror wins
        (instant, offline); otherwise the cache/live API path. Both yield dicts
        with the same keys (trait, gene, beta, se, p, n, tissue, ...)."""
        rows = mirror_lookup(marker)
        if rows is not None:
            return rows  # mirror exists (may be empty list = no associations)
        try:
            payload = self._fetch(marker)
        except Exception:
            return []  # errors surface via status(), not as exceptions here
        fields = payload.get("fields") or []
        return [dict(zip(fields, raw)) for raw in (payload.get("results") or [])]

    def _finding(self, marker: str, row: dict) -> Finding:
        trait = row.get("trait", "unknown trait")
        gene = (row.get("gene") or "").strip() or None
        topic = classify_topic(trait)
        label, kind, accession = humanize_trait(trait)
        # a protein-accession trait reads as a "protein level" measurement
        subject = (f"blood level of protein {label}" if kind == "protein" else label)
        # The sign of the published effect IS the direction of the association, and
        # "linked to age" without it tells the reader nothing about how. Phrased as
        # methylation-relative-to-trait so it stays true whether the study measured a
        # continuous trait or compared two groups. No effect size, no direction claimed.
        beta = row.get("beta")
        try:
            beta = float(beta) if beta is not None and beta != "" else None
        except (TypeError, ValueError):
            beta = None
        if beta:
            way = "lower" if beta < 0 else "higher"
            desc = f"{subject} — associated with {way} methylation at this site"
        else:
            desc = f"linked to {subject}"
        cats = [_TOPIC_CATEGORY.get(topic, Category.TRAIT)]
        detail = {k: row.get(k) for k in
                  ("beta", "se", "p", "n", "tissue", "methylation_array", "chrpos")}
        detail["topic"] = topic
        detail["trait"] = label
        # Study-design variables ("Tissue") describe the sample, not the person.
        # Tagged rather than dropped: the association is a real catalog fact, so
        # it stays in the data and the report decides whether to show it.
        tclass = trait_class(trait)
        if tclass:
            detail["trait_class"] = tclass
        # Resolved from the RAW trait: detail["trait"] below holds the humanized
        # name ("Alpha-2-macroglobulin"), which resolves to no copy at all, so a
        # renderer keying off it would lose every protein trait.
        ckey = trait_copy_key(trait)
        if ckey:
            detail["copy_key"] = ckey
        if gene:
            detail["gene"] = gene
        if accession:
            detail["protein"] = accession
        if row.get("efo"):
            detail["efo"] = row["efo"]     # EFO ontology id (mirror only)
        return Finding(
            marker=marker, source=self.name, description=desc,
            tier=_tier_from(row), categories=cats, detail=detail,
            link=f"https://www.ewascatalog.org/?query={marker}",
            pmids=[str(row["pmid"])] if row.get("pmid") else [])

    def get(self, marker: str) -> list[Finding]:
        return [self._finding(marker, row) for row in self._rows_for(marker)]

    def refresh(self) -> ProviderStatus:
        """Build the local mirror from the EWAS Catalog bulk downloads so lookups
        are instant and offline. Heavy (174 MB download + SQLite build) — meant to
        run on a worker/refresh box, not inline in a request."""
        try:
            summary = build_mirror()
            return ProviderStatus(self.name, Health.OK,
                note=f"mirror built: {summary['n_findings']} associations "
                     f"across {summary['n_studies']} studies")
        except Exception as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE,
                                  note=f"mirror build failed: {e}")

    def status(self) -> ProviderStatus:
        try:
            self._fetch("cg00000029")
            return ProviderStatus(self.name, Health.OK, note="live API reachable")
        except urllib.error.HTTPError as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE, note=f"HTTP {e.code}")
        except Exception as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE,
                                  note=f"{type(e).__name__}: {str(e)[:80]}")
