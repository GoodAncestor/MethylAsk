"""ClinVar provider (NCBI).

Serves the variant-level clinical layer for the SNP (rs) probes carried on
methylation arrays and for probe-overlapping variants (docs/DESIGN.md §3d).

Validated live 2026-07-23 via NCBI E-utilities (allowlisted by default).
Design posture: ClinVar publishes a weekly GRCh38 VCF (~193 MB). The mirror
path (refresh()) downloads and indexes that VCF; the live E-utilities path
here is used for prototyping and single-variant lookups.

Clinical significance is reported as a ClinVar classification, mapped to an
evidence tier — never restated as a diagnosis (see docs/DISCLAIMER.md).
"""
from __future__ import annotations
import json, urllib.parse, urllib.request, urllib.error
from .base import Provider, Finding, Tier, Category, ProviderStatus, Health

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# ClinVar review status -> evidence tier (star rating drives confidence).
_ROBUST = {"practice guideline", "reviewed by expert panel"}
_MODERATE = {"criteria provided, multiple submitters, no conflicts",
             "criteria provided, single submitter"}


def _tier_from_review(status: str) -> Tier:
    s = (status or "").lower()
    if any(k in s for k in _ROBUST):
        return Tier.ROBUST
    if any(k in s for k in _MODERATE):
        return Tier.MODERATE
    return Tier.SPECULATIVE


class ClinVarProvider(Provider):
    name = "clinvar"

    def __init__(self, timeout: int = 30, email: str | None = None):
        self._timeout = timeout
        self._email = email  # supplied via config; NCBI asks for a contact

    def _eutils(self, endpoint: str, params: dict) -> dict:
        if self._email:
            params = {**params, "email": self._email, "tool": "methylask"}
        url = _EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={"User-Agent": "methylask"})
        with urllib.request.urlopen(req, timeout=self._timeout) as r:
            return json.loads(r.read().decode("utf-8", "replace"))

    def get(self, marker: str) -> list[Finding]:
        # marker is an rsID for the SNP probes on the array
        rsid = marker[2:] if marker.lower().startswith("rs") else None
        if not rsid:
            return []
        try:
            es = self._eutils("esearch.fcgi",
                              {"db": "clinvar", "term": f"rs{rsid}", "retmode": "json", "retmax": 5})
            ids = es.get("esearchresult", {}).get("idlist", [])
            if not ids:
                return []
            summ = self._eutils("esummary.fcgi",
                               {"db": "clinvar", "id": ",".join(ids), "retmode": "json"})
        except Exception:
            return []
        out: list[Finding] = []
        result = summ.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            germ = rec.get("germline_classification", {}) or {}
            sig = germ.get("description") or rec.get("clinical_significance", {}).get("description", "unknown")
            review = germ.get("review_status") or ""
            traits = germ.get("trait_set") or []
            trait = traits[0].get("trait_name") if traits else rec.get("title", "variant")
            out.append(Finding(
                marker=marker, source=self.name,
                description=f"ClinVar classifies {marker} as '{sig}' for {trait}",
                tier=_tier_from_review(review),
                categories=[Category.CLINICAL],
                detail={"clinical_significance": sig, "review_status": review},
                link=f"https://www.ncbi.nlm.nih.gov/clinvar/?term=rs{rsid}",
            ))
        return out

    def refresh(self) -> ProviderStatus:
        # TODO: download clinvar.vcf.gz (GRCh38, ~193 MB) -> local index
        return self.status()

    def status(self) -> ProviderStatus:
        try:
            self._eutils("esearch.fcgi",
                        {"db": "clinvar", "term": "rs328", "retmode": "json", "retmax": 1})
            return ProviderStatus(self.name, Health.OK, note="NCBI E-utilities reachable")
        except urllib.error.HTTPError as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE, note=f"HTTP {e.code}")
        except Exception as e:
            return ProviderStatus(self.name, Health.UNAVAILABLE,
                                  note=f"{type(e).__name__}: {str(e)[:80]}")
