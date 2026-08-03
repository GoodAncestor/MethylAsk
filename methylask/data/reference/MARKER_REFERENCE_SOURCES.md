# marker_reference.json — curation notes

Absolute population levels for a marker, quoted from published research. Used to
say "your value is 0.71; never-smokers average 0.85 [PMID]" on a report card.

**These cannot be derived from anything we already hold.** The EWAS Catalog
mirror stores effect sizes — the difference a trait makes to methylation — not
absolute levels in a named reference group. Every entry here is a number read
out of a paper by a person, with the paper recorded.

An unpopulated marker is not a failure: `reference_for()` returns `[]` and the
card reports that no published reference value exists for this marker, which is
a truer statement than a number without a source.

## Schema

```json
{
  "<probe id>": {
    "gene": "<gene symbol, for display only>",
    "references": [
      {
        "group":   "never-smoker",     // the population the value describes
        "beta":    0.85,               // mean beta in that group
        "sd":      0.03,               // OPTIONAL — omit if the paper gives none
        "tissue":  "whole blood",
        "array":   "450k",
        "n":       1793,
        "pmid":    "23691101"
      }
    ]
  }
}
```

`sd` is optional and must be omitted rather than guessed. `describe_position()`
leaves `sigma` as `None` when there is no published SD — a sigma figure derived
from an assumed spread is a fabricated statistic, and these cards are read by
people making health decisions.

Several reference values may exist per marker (e.g. never-smoker, former smoker,
current smoker). List them all; the card can position the sample against each.

## State — first curation pass, 2026-08-03

Populated: 10 rows across 5 probes.

| Probe | Gene | Groups | Source |
|---|---|---|---|
| cg05575921 | AHRR | never / former / current smoker | Zeilinger 2013, PMID 23691101 |
| cg21566642 | 2q37.1 intergenic | never / former / current smoker | Zeilinger 2013, PMID 23691101 |
| cg06690548 | SLC7A11 | AUD cases / healthy controls | Lohoff 2022, PMID 34857913 |
| cg17860381 | NR3C1 | healthy control women | Glad 2017, PMID 28300138 |
| cg03546163 | FKBP5 | birth-cohort infants | Mulder 2017, PMID 28401840 |

Carry these caveats into any copy that quotes them:

- **The smoking values are medians, not means** (`stat: "median"`; the paper
  prints only median + IQR). Do not call them averages and do not derive a sigma.
- **cg17860381 is n=16**, a clinical study's control arm — not a population.
- **cg03546163 is cord blood**, i.e. neonatal. Comparing an adult sample to it
  requires saying so.
- **cg06690548 is alcohol-use-disorder cases vs controls**, a heavy-vs-low proxy.
  No paper publishes per-intake-band means, so a claim like "vs light drinkers"
  is not supported.
- Only one row (cg03546163) has a published SD. Assume `sd` is absent.

## Probe IDs that do not exist

Checked against the bundled Zhou-lab manifests via `normalize.Manifest` — these
are on **no** array generation (27k/450k/EPIC/EPICv2) and are not curatable:

    cg20241083    cg09649689    cg11761638

They came in with the Cellscript prototype's marker list. `cg16232979` *is* a
real probe but maps to **chr19:16076820**, which is not CLOCK (chr4), PER2
(chr2) or ARNTL (chr11) — so its attribution to a circadian gene is also wrong.
Validate any new marker against the manifests before commissioning research on it.

## Composites — do not exist as specified

Neither the "19-site PM2.5" nor the "7-site early-life adversity" composite
exists in the literature as described. The only 19-CpG air-pollution set is a
7 + 12 split that is **NO2**-associated, never summed into a published score.
The well-known "7-CpG" signature is the Vidal-Bralo/BASE-II **chronological-age
clock** — a name collision, unrelated to adversity. Real adversity composites
have 3, 4, 9 or 14 CpGs and none publishes a score mean + SD.

Adopting a composite is a product decision, not a research gap. Any composite we
adopt needs its published probe list **and weights**, not just a distribution —
a distribution for a score computed differently is meaningless.
