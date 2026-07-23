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

## Layout

    methylask/
      providers/    data-layer: one provider per reference database
      ingest/       file-format parsers (CSV/GEO, IDAT, bedMethyl, ...)
      report/       report model + HTML/PDF rendering
    data/reference/ small static reference files committed to the repo (<100 MB)
    scripts/        refresh + build CLIs
    docs/           design, validation, disclaimer

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
