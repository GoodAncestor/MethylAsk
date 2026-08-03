"""Replication context for a marker, derived from EWAS mirror rows.

The robustness statement on a finding is computed from the evidence itself —
how many independent studies saw it, in how many people, in which tissues, and
whether they agreed on the direction of effect — rather than assigned by hand.
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class ReplicationContext:
    direction: str          # "consistent" | "mixed"
    n_studies: int          # distinct publications, not rows
    n_participants: int     # people, counted once per study
    tissues: list[str]      # tissues the association was actually studied in
    tissue_supported: bool  # sample's tissue appears in that literature


# A marker studied in "whole blood" supports a blood sample; the mirror's tissue
# strings are free text, so matching is on the containing word, not equality.
def tissue_matches(sample: str, studied: list[str]) -> bool:
    s = sample.strip().lower()
    return any(s in t.lower() or t.lower() in s for t in studied)


_tissue_match = tissue_matches   # internal alias, kept for readability below


def summarize_replication(rows: list[dict],
                          sample_tissue: str | None = None) -> ReplicationContext:
    """Summarise per-association rows (as returned by ewas_mirror.mirror_lookup).

    A single publication contributes one row per trait it reported, so studies
    and participants are collapsed by pmid — summing the n column row-wise would
    inflate a well-replicated marker by however many traits each paper measured.
    """
    signs = {(1 if r["beta"] > 0 else -1) for r in rows if r.get("beta") is not None}

    per_study: dict[str, int] = {}
    for r in rows:
        pmid = r.get("pmid") or ""
        n = r.get("n") or 0
        per_study[pmid] = max(per_study.get(pmid, 0), int(n))

    tissues = sorted({(r.get("tissue") or "").strip()
                      for r in rows if (r.get("tissue") or "").strip()})

    return ReplicationContext(
        direction="mixed" if len(signs) > 1 else "consistent",
        n_studies=len(per_study),
        n_participants=sum(per_study.values()),
        tissues=tissues,
        tissue_supported=(_tissue_match(sample_tissue, tissues)
                          if sample_tissue else True),
    )
