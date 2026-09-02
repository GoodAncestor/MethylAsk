# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Fold a marker's EWAS rows into one finding per trait.

The catalogue has one row per publication, marker, and trait. An aggregated
finding carries the study count, participant count, tissues, and direction.
"""
from __future__ import annotations

from biocore.providers.base import Finding, Tier

from .evidence import summarize_replication


_ORDER = {
    Tier.ROBUST: 0,
    Tier.MODERATE: 1,
    Tier.SPECULATIVE: 2,
    Tier.UNKNOWN: 3,
}


def _number(value):
    try:
        return None if value in (None, "") else float(value)
    except (TypeError, ValueError):
        return None


def aggregate_by_trait(
    findings: list[Finding], sample_tissue: str | None
) -> list[Finding]:
    """Collapse EWAS rows by marker and trait while preserving other findings."""
    groups: dict[tuple[str, str], list[Finding]] = {}
    order: list[tuple[str, str]] = []
    out: list[Finding] = []

    for finding in findings:
        if finding.source != "ewas_catalog":
            out.append(finding)
            continue
        key = (finding.marker, str((finding.detail or {}).get("trait") or ""))
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(finding)

    for key in order:
        grouped = groups[key]
        if len(grouped) == 1:
            out.append(grouped[0])
            continue

        rows = [
            {
                "pmid": (finding.pmids or [""])[0],
                "beta": _number(finding.detail.get("beta")),
                "p": _number(finding.detail.get("p")),
                "n": int(_number(finding.detail.get("n")) or 0),
                "tissue": finding.detail.get("tissue"),
            }
            for finding in grouped
        ]
        context = summarize_replication(rows, sample_tissue)
        best = min(
            grouped,
            key=lambda finding: (
                _ORDER[finding.tier],
                _number(finding.detail.get("p")) or 1.0,
            ),
        )
        detail = dict(best.detail)
        detail.update(
            {
                "n_studies": context.n_studies,
                "n_participants": context.n_participants,
                "direction": context.direction,
                "tissues": context.tissues,
                "tissue_supported": (
                    context.tissue_supported if sample_tissue else None
                ),
                "rows": rows,
                "p": min(
                    (row["p"] for row in rows if row["p"] is not None),
                    default=None,
                ),
                "n": context.n_participants,
            }
        )
        pmids = sorted(
            {pmid for finding in grouped for pmid in (finding.pmids or [])}
        )
        out.append(
            Finding(
                marker=best.marker,
                source="ewas_catalog",
                description=best.description,
                tier=best.tier,
                categories=list(best.categories),
                detail=detail,
                link=best.link,
                pmids=pmids,
            )
        )
    return out
