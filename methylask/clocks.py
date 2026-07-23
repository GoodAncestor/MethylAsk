"""Epigenetic clocks — run locally from published coefficients.

A clock is a weighted sum over a fixed CpG set (docs/DESIGN.md §4.3, §8.1):
each clock ships as a Probe/Coefficient table in data/reference/clocks/.
Computation is a dot product; nothing leaves the server.

Two age transforms are needed (see clocks/SOURCES.md):
  - linear:            age = intercept + sum(coef_i * beta_i)      (Hannum, PhenoAge)
  - anti_log_linear:   Horvath's transform on the linear predictor  (Horvath 2013 / Skin&Blood)

Coverage matters: some clock CpGs are absent from EPIC/EPICv2 (docs/DESIGN.md §3c).
The engine reports how many CpGs were found so a low-coverage estimate is flagged
rather than silently biased.
"""
from __future__ import annotations
import csv, math
from dataclasses import dataclass
from pathlib import Path

_CLOCK_DIR = Path(__file__).parent / "data" / "reference" / "clocks"

# transform + adult-age constant per clock file (matches clocks/SOURCES.md)
_CLOCKS = {
    "Horvath2013_PanTissue": ("anti_log_linear", 20.0),
    "Horvath2018_SkinBlood": ("anti_log_linear", 20.0),
    "Hannum2013_Blood":      ("linear",          None),
    "Levine2018_PhenoAge":   ("linear",          None),
}


def _anti_trafo(x: float, adult_age: float) -> float:
    """Inverse of Horvath's age transformation (Horvath 2013, Genome Biology)."""
    if x < 0:
        return (1.0 + adult_age) * math.exp(x) - 1.0
    return (1.0 + adult_age) * x + adult_age


@dataclass
class ClockResult:
    clock: str
    age: float | None            # predicted DNAm age (years), None if unusable
    n_cpg: int                   # CpGs in the clock
    n_found: int                 # CpGs present in the sample
    coverage: float              # n_found / n_cpg
    low_coverage: bool           # True if coverage below the flag threshold

    @property
    def note(self) -> str:
        pct = f"{self.coverage*100:.0f}%"
        if self.low_coverage:
            return f"{self.n_found}/{self.n_cpg} CpGs ({pct}) — low coverage, estimate may be biased"
        return f"{self.n_found}/{self.n_cpg} CpGs ({pct})"


class Clock:
    def __init__(self, name: str, min_coverage: float = 0.80):
        if name not in _CLOCKS:
            raise ValueError(f"unknown clock {name!r}; have {list(_CLOCKS)}")
        self.name = name
        self.transform, self.adult_age = _CLOCKS[name]
        self.min_coverage = min_coverage
        self.intercept = 0.0
        self.weights: dict[str, float] = {}
        self._load()

    def _load(self) -> None:
        with open(_CLOCK_DIR / f"{self.name}.csv", newline="") as fh:
            for row in csv.DictReader(fh):
                probe, coef = row["Probe"], float(row["Coefficient"])
                if probe == "Intercept":
                    self.intercept = coef
                else:
                    self.weights[probe] = coef

    @property
    def n_cpg(self) -> int:
        return len(self.weights)

    def predict(self, betas: dict[str, float]) -> ClockResult:
        """betas keyed by BASE probe id (suffix-stripped, see normalize.base_probe)."""
        acc, found = self.intercept, 0
        for probe, w in self.weights.items():
            b = betas.get(probe)
            if b is not None:
                acc += w * b
                found += 1
        coverage = found / self.n_cpg if self.n_cpg else 0.0
        low = coverage < self.min_coverage
        if found == 0:
            return ClockResult(self.name, None, self.n_cpg, 0, 0.0, True)
        age = _anti_trafo(acc, self.adult_age) if self.transform == "anti_log_linear" else acc
        return ClockResult(self.name, age, self.n_cpg, found, coverage, low)


def available() -> list[str]:
    return list(_CLOCKS)


def run_all(betas: dict[str, float], min_coverage: float = 0.80) -> list[ClockResult]:
    return [Clock(n, min_coverage).predict(betas) for n in _CLOCKS]
