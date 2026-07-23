# MethylAsk

A web tool that reads a person's methylation/epigenomic file, matches the markers
against public scientific databases, and produces a readable report of what is
known about them — clinical relevance and popular-interest findings side by side,
each tagged with an evidence tier.

Runs on a local server or private cloud so personal data stays under the
operator's control.

## Status

Prototype scaffold. See `docs/DESIGN.md` for the full architecture and
`docs/VALIDATION.md` for live data-source validation results.

## Built on bio-core

MethylAsk depends on [bio-core](https://github.com/GoodAncestor/bio-core) for
all organism-agnostic mechanism — the provider interface and registry, evidence
tiering, the HTML/PDF report renderer, and the resumable fetch helper. MethylAsk
itself holds only the human-methylation **knowledge**: the methylation databases
it queries, the epigenetic clocks, and the Illumina-array normalization.

    methylask/
      providers/    methylation-specific providers (EWAS Catalog, ClinVar, GDC)
                    — each imports the shared Provider/Finding/Tier from biocore
      ingest/       file-format parsers (CSV/GEO, IDAT, bedMethyl, ...)
      normalize.py  probe/rsID → canonical genome coordinate (bundled manifests)
      clocks.py     epigenetic-clock engine (Horvath, Hannum, PhenoAge, ...)
      cli.py        `methylask status|refresh|report`
    data/reference/ small static reference files committed to the repo (<100 MB)
    scripts/        refresh + build CLIs
    docs/           design, validation, disclaimer

The provider registry, evidence tiering, and report renderer live in bio-core
(`biocore.providers`, `biocore.report`) and are imported, not duplicated.

## Quick start

    pip install -e .
    methylask status                     # provider health + cache ages
    methylask refresh --provider all     # build/refresh local caches
    methylask report sample.csv --pdf    # produce a report from a sample

## Data

Small static reference files (array manifests, clock coefficients) live in this
repo. Large corpora (GDC/TCGA methylation, ClinVar, EWAS Catalog dump) are
mirrored to local server disk by `methylask refresh` — see `docs/DESIGN.md` §3.4.

## Disclaimer

MethylAsk reports research associations, not medical diagnoses. See
[docs/DISCLAIMER.md](docs/DISCLAIMER.md) — the single source of truth for all
disclaimer language.
