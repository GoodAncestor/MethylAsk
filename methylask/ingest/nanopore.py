"""Map native ONT methylation to array marker IDs without array calibration.

Coordinates are always 0-based forward-reference CpG starts. The hg38 Zhou
manifests supply genomic addresses, not evidence that an ONT read fraction is
interchangeable with an Illumina beta. No imputation, liftover, clocks, or
reference-population classification is performed here.
"""
from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass, field
from typing import Iterable

from biocore.io.bedmethyl import read_sites
from biocore.io.modbam import pileup_methyl_targets
from biocore.methylation.model import Context, MethylSite

from .beta_matrix import Sample
from ..normalize import _FILES, _MANIFEST_DIR, base_probe

PLATFORM_WARNING = (
    "ONT fractions are native sequencing measurements, not calibrated Illumina betas. "
    "Array-trained clocks and reference-group classification are not validated for this "
    "conversion and must remain disabled unless explicitly run as exploratory analyses. "
    "5mC and 5hmC are kept separate; bisulfite-array values can include both."
)
_PRIMARY_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}


@dataclass
class NanoporeSample(Sample):
    coverage: dict[str, int] = field(default_factory=dict)
    stats: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=lambda: [PLATFORM_WARNING])


def _chrom(chrom: str) -> str:
    c = chrom.removeprefix("chr")
    return "chr" + c if c in _PRIMARY_CHROMS else chrom


def _probe_map(array: str, supplied: dict[str, tuple[str, int]] | None):
    if supplied is not None:
        source = "caller-supplied probe_map (0-based forward CpG start)"
        items = supplied.items()
    else:
        if array not in _FILES:
            raise ValueError(f"Unknown array manifest: {array}")
        path = _MANIFEST_DIR / _FILES[array]
        source = path.name
        def rows():
            with gzip.open(path, "rt") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if (not row.get("Probe_ID", "").startswith("cg")
                            or row.get("CpG_chrm", "NA") == "NA"
                            or row.get("CpG_beg", "NA") == "NA"):
                        continue
                    yield row["Probe_ID"], (row["CpG_chrm"], int(row["CpG_beg"]))
        items = rows()
    mapping, ambiguous = {}, set()
    for probe, (chrom, pos) in items:
        probe = base_probe(probe)
        if not probe.startswith("cg") or not chrom or not isinstance(pos, int) or pos < 0:
            raise ValueError("Probe map must contain CpG IDs and nonnegative 0-based coordinates")
        locus = (_chrom(chrom), pos)
        if probe in mapping and mapping[probe] != locus:
            ambiguous.add(probe)
        mapping[probe] = locus
    for probe in ambiguous:
        mapping.pop(probe)
    if not mapping:
        raise ValueError("No unambiguous CpG coordinates in the selected manifest")
    reverse = {}
    for probe, locus in mapping.items():
        reverse.setdefault(locus, []).append(probe)
    return reverse, source, len(ambiguous)


def _prepare(reference_build, array, min_coverage, mod_code, probe_map):
    if reference_build not in {"hg38", "GRCh38"}:
        raise ValueError("ONT marker mapping requires explicit hg38/GRCh38; no automatic liftover")
    if not isinstance(min_coverage, int) or min_coverage < 1:
        raise ValueError("min_coverage must be a positive integer")
    if mod_code not in {"m", "h"}:
        raise ValueError("Choose mod_code='m' (5mC) or 'h' (5hmC); combined codes are not accepted")
    return _probe_map(array, probe_map)


def _convert(sites: Iterable[MethylSite], *, reverse, source, reference_build,
             array, min_coverage, mod_code, combined_strands, format_name,
             stats, provenance, ambiguous):
    sample = NanoporeSample(meta={"source": source, "platform": "ONT", "array": array,
                                  "build": "hg38", "modification": mod_code,
                                  "measurement": "native_modified_read_fraction"})
    tally, seen = {}, set()
    for site in sites:
        if site.mod_code != mod_code or site.context not in {Context.CG, Context.UNKNOWN}:
            continue
        if site.strand == "." and not combined_strands:
            raise ValueError("Unstranded bedMethyl requires combined_strands=True and forward CpG-start coordinates")
        if combined_strands and site.strand != ".":
            raise ValueError("combined_strands=True requires exclusively '.' strand rows")
        pos = site.pos - 1 if site.strand == "-" else site.pos
        locus = (_chrom(site.chrom), pos)
        if locus not in reverse:
            continue
        key = (*locus, site.strand)
        if key in seen:
            raise ValueError(f"Duplicate {mod_code} measurement at {locus[0]}:{pos} strand {site.strand}; "
                             "do not concatenate overlapping motif, haplotype, or replicate pileups")
        seen.add(key)
        cell = tally.setdefault(locus, [0, 0, 0])
        cell[0] += site.n_mod
        cell[1] += site.coverage
        cell[2] += site.n_other_mod
    observed = low = 0
    for locus, (nmod, cov, _) in tally.items():
        for probe in reverse[locus]:
            sample.coverage[probe] = cov
            observed += 1
            if cov >= min_coverage:
                sample.betas[probe] = nmod / cov
            else:
                low += 1
    total = sum(len(probes) for probes in reverse.values())
    sample.stats = {**stats, "total_loci": len(reverse), "observed_loci": len(tally),
                    "qualifying_loci": sum(cov >= min_coverage for _, cov, _ in tally.values()),
                    "total_probes": total, "observed_probes": observed,
                    "mapped_probes": len(sample.betas), "low_coverage_probes": low,
                    "missing_probes": total - observed, "ambiguous_manifest_probes": ambiguous,
                    "min_coverage": min_coverage,
                    "status": "ready" if sample.betas else "no_qualifying_markers"}
    sample.provenance = {**provenance, "format": format_name, "source": source,
                         "reference_build": "GRCh38", "declared_reference_build": reference_build,
                         "coordinate_system": "0-based forward CpG start",
                         "mod_code": mod_code, "fraction_denominator": "all valid calls including other modifications",
                         "strand_policy": "already combined" if combined_strands else "sum opposite CpG cytosines",
                         "min_coverage": min_coverage, "imputation": "none",
                         "platform_validation": "unvalidated_cross_platform",
                         "array_clocks_allowed_by_default": False}
    if not sample.betas:
        sample.warnings.append("No mapped CpG passed the coverage threshold; absent/low-coverage loci are not zero methylation.")
    if ambiguous:
        sample.warnings.append(f"Excluded {ambiguous} probe IDs with conflicting manifest coordinates.")
    return sample


def read_bedmethyl(path: str, *, reference_build: str, array: str = "EPICv2",
                  min_coverage: int = 5, mod_code: str = "m",
                  combined_strands: bool = False,
                  probe_map: dict[str, tuple[str, int]] | None = None) -> NanoporeSample:
    """Stream modkit 18-column bedMethyl, retaining only mapped target loci.

    Standard single-letter modification names and motif-annotated names are
    accepted. The declared build is required because BED has no assembly header.
    Duplicate selected-modification rows at target sites fail instead of inflating
    coverage. Opposite strands sum counts before applying min_coverage.
    """
    reverse, manifest, ambiguous = _prepare(reference_build, array, min_coverage, mod_code, probe_map)
    stats = {}
    sites = read_sites(path, strict=True, mod_codes={mod_code}, stats=stats)
    return _convert(sites, reverse=reverse, source=str(path), reference_build=reference_build,
                    array=array, min_coverage=min_coverage, mod_code=mod_code,
                    combined_strands=combined_strands, format_name="bedMethyl", stats=stats,
                    provenance={"manifest": manifest, "sample_identity": "not encoded in bedMethyl"},
                    ambiguous=ambiguous)


def read_modbam(path: str, *, reference_build: str, reference_fasta: str,
                array: str = "EPICv2", min_coverage: int = 5, mod_code: str = "m",
                min_prob: float = 0.8, min_mapping_quality: int = 20,
                probe_map: dict[str, tuple[str, int]] | None = None) -> NanoporeSample:
    """Read a single-sample aligned modBAM, retaining only manifest targets.

    The reference must be the declared GRCh38 alignment FASTA with an index.
    Unsorted BAM is supported via a single sequential pass. Missing/ambiguous
    MM/ML calls and filtered alignments do not contribute valid coverage.
    """
    import pysam
    reverse, manifest, ambiguous = _prepare(reference_build, array, min_coverage, mod_code, probe_map)
    stats, provenance = {}, {"manifest": manifest}
    with pysam.AlignmentFile(path, "rb") as bam:
        contigs = {_chrom(chrom): chrom for chrom in bam.references}
        if len(contigs) != len(bam.references):
            raise ValueError("BAM contains ambiguous chromosome aliases")
    targets = {}
    for chrom, pos in reverse:
        if chrom in contigs:
            targets.setdefault(contigs[chrom], set()).update((pos, pos + 1))
    sites = pileup_methyl_targets(path, targets=targets, reference_fasta=reference_fasta,
                                  mod_code=mod_code, min_prob=min_prob,
                                  min_mapping_quality=min_mapping_quality,
                                  stats=stats, provenance=provenance)
    return _convert(sites, reverse=reverse, source=str(path), reference_build=reference_build,
                    array=array, min_coverage=min_coverage, mod_code=mod_code, combined_strands=False,
                    format_name="modBAM", stats=stats, provenance=provenance, ambiguous=ambiguous)
