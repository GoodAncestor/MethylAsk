# Native Oxford Nanopore methylation bridge

`methylask.ingest.nanopore` converts single-sample aligned modBAM or modkit
18-column bedMethyl into a `NanoporeSample`. This extends the existing `Sample`
interface (`betas`, `markers`, `meta`) with `coverage`, `stats`, `provenance`,
and `warnings`. The `betas` field is an API compatibility name: these values
are **native modified-read fractions, not calibrated Illumina beta values**.

```python
from methylask.ingest.nanopore import read_bedmethyl, read_modbam

sample = read_bedmethyl(
    "sample.cpg.bed.gz", reference_build="GRCh38", min_coverage=5,
    mod_code="m", combined_strands=False,
)
# Alternatively, read aligned MM/ML calls directly, using its indexed FASTA:
sample = read_modbam(
    "sample.mod.bam", reference_build="GRCh38",
    reference_fasta="GRCh38.fa", min_coverage=5, min_prob=0.8,
)
print(sample.stats, sample.provenance)
```

## Contracts

- Build is mandatory and must be `hg38` or `GRCh38`. No inferred build or
  liftover. BedMethyl cannot verify its own assembly: the supplied declaration
  must come from the alignment provenance. BAM/FASTA contig lengths are checked,
  but matching lengths alone do not authenticate assembly identity.
- The bundled Zhou EPICv2 hg38 manifest is the default coordinate map; HM450
  and EPIC are selectable with `array=`. Coordinates use `CpG_beg`, the 0-based
  forward-reference CpG start. Replicate suffixes collapse; conflicting probe
  coordinates are excluded. `probe_map={"cg...": ("chr1", 10)}` supports
  caller-supplied coordinates and hermetic tests with the same convention.
- A plus-strand C at *p* and minus-strand C at *p+1* describe the same CpG;
  their counts sum before the coverage threshold. `chr1`/`1` are aliases.
  Unstranded rows require explicit `combined_strands=True` and must already
  have forward CpG-start coordinates. Mixed combined/separate strand rows and
  duplicate selected-modification target rows fail rather than inflate coverage.
- `mod_code="m"` is 5mC; `"h"` is 5hmC. Values are selected `N_mod / N_valid`.
  Other confident modifications remain in the denominator. The two channels
  are never added implicitly. Plain modification names and motif names such as
  `m,CG,0` are accepted. Non-CpG contexts are excluded.
- BedMethyl is streamed once, with only target loci retained. Invalid row
  structure/count totals fail with a line number. Zero-coverage and absent loci
  are never interpreted as zero methylation or imputed. Coverage for observed
  probes includes those below the threshold; `stats` separates missing from
  low coverage and reports mapped probes, target/observed/qualifying loci.
- Direct BAM mode retains only manifest target cytosines. It accepts unsorted
  BAM in a single pass; memory scales with target loci and the longest read.
  It excludes unmapped, secondary, supplementary, duplicate, QC-failed, and
  low-MAPQ alignments. Reference sequence verifies CpG context. Multiple BAM
  `RG:SM` identities fail; an absent sample ID is recorded as an empty list.
- Direct BAM mode selects the most probable state among canonical and all
  modeled cytosine modifications, using ML bin midpoints. The winning state
  must exceed the specified confidence threshold. Unknown ML values and
  `?`-omitted calls do not count; `.`/default omissions retain canonical calls.
  This is a transparent research estimator, **not a claim of modkit threshold
  equivalence**. Use modkit-produced bedMethyl for the main calling workflow
  and preserve its caller versions, models, filters and sample identity in the
  enclosing run manifest.

## Interpretation boundary

Array-trained clocks and nearest reference-group classification must be
disabled by default for these inputs. Conversion does not validate those
models, correct tissue effects, or estimate clinical risk. Bisulfite arrays
can measure 5mC and 5hmC together, whereas this API keeps them separate.
Exploratory model use requires explicit selection and prominent qualification.
EWAS annotations describe published associations; they are not personal
diagnoses or validated comparisons of the ONT fraction against an array cohort.

Whole-profile reference annotation should use local mirrors: constructors for
`EwasCatalogProvider`, `GdcProvider`, and `ClinVarProvider` accept
`offline_only=True`. Status and lookups then cannot trigger live fallback.
EWAS bulk lookup uses one read-only SQLite connection and bounded SQL batches.
Unavailable mirrors produce explicit provider status, not a long series of
web requests. The methylation ClinVar provider currently has no local backend;
variant interpretation belongs to the separate variant pipeline.

## Sources and validation

- [ONT modkit bedMethyl contract](https://github.com/nanoporetech/modkit#description-of-bedmethyl-output).
- [SAM MM/ML specification](https://github.com/samtools/hts-specs/blob/master/SAMtags.tex).
- [pysam modified base and reference coordinate API](https://pysam.readthedocs.io/en/latest/api.html).

Tests use generated sequences/calls and tiny marker maps, including opposite
strands, mixed 5mC/5hmC states, unknown/ambiguous calls, canonical omissions,
filtered/duplicate reads, count corruption, missingness and cross-format count
agreement. They establish conversion behavior only; no paired ONT/array samples
or clinical validation have been run.
