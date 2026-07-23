"""Beta-value matrix ingest (the ship-first format, docs/DESIGN.md §4.1).

Reads a CSV/TSV where rows are CpG probe ids (cg########) and one column holds
the beta value for the sample. Produces the common internal object: a mapping
of marker id -> beta, plus a metadata header.

This is a *format conversion*, not normalization — the betas are already
processed. Raw IDAT -> beta is a separate, heavier step (docs/DESIGN.md §8.3).
"""
from __future__ import annotations
import csv
from dataclasses import dataclass, field


@dataclass
class Sample:
    betas: dict[str, float] = field(default_factory=dict)  # probe id -> beta
    meta: dict[str, str] = field(default_factory=dict)     # array, build, source

    @property
    def markers(self) -> list[str]:
        return list(self.betas)


def read_beta_matrix(path: str, probe_col: int = 0, value_col: int = 1,
                     delimiter: str | None = None) -> Sample:
    with open(path, newline="") as fh:
        if delimiter is None:
            head = fh.readline()
            delimiter = "\t" if head.count("\t") >= head.count(",") else ","
            fh.seek(0)
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader, None)
        s = Sample(meta={"source": path})
        for row in reader:
            if len(row) <= max(probe_col, value_col):
                continue
            probe = row[probe_col].strip()
            if not probe.startswith(("cg", "ch", "rs")):
                continue
            try:
                s.betas[probe] = float(row[value_col])
            except ValueError:
                continue
        return s
