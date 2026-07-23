"""Normalization — resolve every marker to a canonical (probe, chrom:pos, build).

The stage that prevents confidently-wrong answers (docs/DESIGN.md §4.2). The
same physical CpG is named differently across array generations and sits at
different coordinates across genome builds. Before any database lookup we map
each probe to a canonical location using the bundled Zhou-lab manifests.

EPICv2 wrinkle: probe IDs carry replicate suffixes (cg#####_TC21). The base id
(before the underscore) is what maps to a location and what older-array data is
keyed by, so we split it off and keep both.

A probe that cannot be confidently mapped is FLAGGED and excluded from
interpretation rather than guessed at.
"""
from __future__ import annotations
import gzip, csv
from dataclasses import dataclass, field
from pathlib import Path

_MANIFEST_DIR = Path(__file__).parent / "data" / "reference" / "manifests"
_FILES = {
    "HM450": "HM450.hg38.manifest.tsv.gz",
    "EPIC": "EPIC.hg38.manifest.tsv.gz",
    "EPICv2": "EPICv2.hg38.manifest.tsv.gz",
}


def base_probe(probe: str) -> str:
    """Strip EPICv2 replicate suffix: 'cg00000029_TC21' -> 'cg00000029'."""
    return probe.split("_", 1)[0]


@dataclass
class CanonicalMarker:
    probe: str                 # original id as supplied
    base: str                  # base probe id (suffix stripped)
    chrom: str | None          # e.g. 'chr16'
    pos: int | None            # CpG start, hg38
    build: str = "hg38"
    mapped: bool = True

    @property
    def locus(self) -> str | None:
        return f"{self.chrom}:{self.pos}" if self.mapped else None


class Manifest:
    """Bundled probe -> hg38 coordinate map for one array generation."""

    def __init__(self, array: str):
        if array not in _FILES:
            raise ValueError(f"unknown array {array!r}; expected one of {list(_FILES)}")
        self.array = array
        self._map: dict[str, tuple[str, int]] = {}
        self._load()

    def _load(self) -> None:
        path = _MANIFEST_DIR / _FILES[self.array]
        with gzip.open(path, "rt") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                pid = row.get("Probe_ID") or ""
                chrom = row.get("CpG_chrm")
                beg = row.get("CpG_beg")
                if not pid or not chrom or chrom == "NA" or not beg or beg == "NA":
                    continue
                self._map[base_probe(pid)] = (chrom, int(beg))

    def __len__(self) -> int:
        return len(self._map)

    def resolve(self, probe: str) -> CanonicalMarker:
        b = base_probe(probe)
        hit = self._map.get(b)
        if hit is None:
            return CanonicalMarker(probe, b, None, None, mapped=False)
        return CanonicalMarker(probe, b, hit[0], hit[1])


@dataclass
class NormalizationResult:
    mapped: dict[str, CanonicalMarker] = field(default_factory=dict)
    unmapped: list[str] = field(default_factory=list)

    def summary(self) -> str:
        total = len(self.mapped) + len(self.unmapped)
        return f"{len(self.mapped)}/{total} probes mapped, {len(self.unmapped)} flagged unmapped"


def normalize(probes: list[str], array: str) -> NormalizationResult:
    mani = Manifest(array)
    res = NormalizationResult()
    for p in probes:
        cm = mani.resolve(p)
        if cm.mapped:
            res.mapped[p] = cm
        else:
            res.unmapped.append(p)
    return res
