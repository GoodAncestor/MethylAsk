"""Published reference values for a marker, quoted with their source.

A card that says "0.71 beta, never-smokers average 0.85" is making a claim about
a population, and that claim belongs to a paper. This table holds those quoted
values keyed by probe, each with the PMID it came from, so the report cites
rather than asserts.

Why this is curated and not computed: the EWAS Catalog mirror stores effect
sizes — the difference a trait makes to methylation — not absolute levels in a
reference group. There is no arithmetic that turns one into the other. A marker
we surface without a published reference must therefore report the absence,
which is why reference_for() returns an empty list instead of raising.
"""
from __future__ import annotations
import json
from dataclasses import dataclass
from pathlib import Path

DEFAULT_TABLE = Path(__file__).parent / "data" / "reference" / "marker_reference.json"


@dataclass
class Position:
    """Where a sample sits relative to one published reference value."""
    reference_group: str
    reference_beta: float
    delta: float                 # sample - reference, in beta units
    pmid: str
    sigma: float | None = None   # only when the source published an SD


def load_reference_table(path: Path | None = None) -> dict:
    path = Path(path or DEFAULT_TABLE)
    if not path.exists():
        return {}
    with open(path) as fh:
        return json.load(fh)


def reference_for(table: dict, marker: str) -> list[dict]:
    """Published reference values for a marker; [] when none is curated."""
    return list(table.get(marker, {}).get("references", []))


def describe_position(sample_beta: float, ref: dict) -> Position:
    """Position a sample against one reference value.

    sigma is left None unless the source published an SD — deriving a
    sigma figure without one would be presenting a fabricated statistic.
    """
    ref_beta = float(ref["beta"])
    sd = ref.get("sd")
    return Position(
        reference_group=ref["group"],
        reference_beta=ref_beta,
        delta=sample_beta - ref_beta,
        pmid=str(ref.get("pmid", "")),
        sigma=((sample_beta - ref_beta) / float(sd)) if sd else None,
    )
