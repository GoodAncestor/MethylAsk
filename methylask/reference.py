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
    stat: str | None = None      # "median" / "mean" / "adjusted mean", as published
    tissue: str | None = None    # cord blood and buccal are not adult whole blood
    n: int | None = None         # an n of 16 is not a population reference


def load_reference_table(path: Path | None = None) -> dict:
    """Load and validate the curated table.

    Betas are proportions. A value outside 0-1 means a percentage (85) or an
    effect size (-0.17) was pasted into a beta field — a silent, large error, so
    it fails the load rather than reaching a report.
    """
    path = Path(path or DEFAULT_TABLE)
    if not path.exists():
        return {}
    with open(path) as fh:
        table = json.load(fh)
    for marker, entry in table.items():
        for ref in entry.get("references", []):
            beta = float(ref["beta"])
            if not 0.0 <= beta <= 1.0:
                raise ValueError(
                    f"{marker}: beta {beta} is outside 0-1 — a beta is a "
                    f"proportion, so this looks like a percentage or an effect size")
    return table


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
        stat=ref.get("stat"),
        tissue=ref.get("tissue"),
        n=ref.get("n"),
    )
