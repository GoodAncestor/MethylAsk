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


_MEANING_FIELDS = ("label", "what_was_read", "what_it_is_not")


def marker_meaning(table: dict, marker: str) -> dict:
    """Plain-language copy for a marker: what it is, and what it is not.

    A probe id and a beta value mean nothing to the person reading the report.
    Each curated marker carries a human label plus the pair of statements that
    make a finding usable — what the measurement is, and the limits of what it
    can support. Empty dict when the marker is not curated.
    """
    entry = table.get(marker, {})
    return {k: entry[k] for k in _MEANING_FIELDS if entry.get(k)}


def reference_for(table: dict, marker: str) -> list[dict]:
    """Published reference values for a marker; [] when none is curated."""
    return list(table.get(marker, {}).get("references", []))


@dataclass
class MarkerPositions:
    """Every published reference for one probe, positioned — or withheld."""
    probe: str
    sample_beta: float
    positions: list[Position]
    suppressed_reason: str | None = None


def positions_for_sample(betas: dict[str, float], tissue: str | None = None,
                         table: dict | None = None) -> list[MarkerPositions]:
    """Position a whole sample against the curated table.

    Local arithmetic — no network — so this runs on the full beta profile rather
    than a capped subset, the same way the clocks do.

    A reference measured in a tissue the sample is not is WITHHELD, not shown
    with a caveat. The pediatric buccal demo is why: at cg05575921 it scores
    0.56, whose nearest whole-blood neighbour is the current-smoker median, so a
    child's cheek swab would read as a smoker. This mirrors clocks.py, where a
    blood-trained clock on buccal is marked not-valid rather than displayed.
    """
    from .evidence import tissue_matches
    table = load_reference_table() if table is None else table
    out: list[MarkerPositions] = []
    for probe, beta in betas.items():
        refs = reference_for(table, probe)
        if not refs:
            continue
        ref_tissues = sorted({r.get("tissue", "") for r in refs if r.get("tissue")})
        if tissue and ref_tissues and not tissue_matches(tissue, ref_tissues):
            out.append(MarkerPositions(
                probe, beta, [],
                suppressed_reason=(
                    f"this sample is {tissue}; the only published reference values "
                    f"for {probe} are {', '.join(ref_tissues)}, and methylation "
                    f"levels are not comparable across tissues")))
            continue
        out.append(MarkerPositions(
            probe, beta, [describe_position(beta, r) for r in refs]))
    return out


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
