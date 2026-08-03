# marker_reference.json — curation notes

Absolute population levels for a marker, quoted from published research. Used to
say "your value is 0.71; never-smokers average 0.85 [PMID]" on a report card.

**These cannot be derived from anything we already hold.** The EWAS Catalog
mirror stores effect sizes — the difference a trait makes to methylation — not
absolute levels in a named reference group. Every entry here is a number read
out of a paper by a person, with the paper recorded.

The table ships empty. An unpopulated marker is not a failure: `reference_for()`
returns `[]` and the card reports that no published reference value exists for
this marker, which is a truer statement than a number without a source.

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

## Wanted

Markers surfaced by the Cellscript prototype, each needing a sourced reference
value before its card can quote one. Morgan's prototype already displays numbers
for the starred ones, so her sources are the fastest path.

| Marker | Gene / signature | Value needed | Prototype shows |
|---|---|---|---|
| cg05575921 ★ | AHRR — tobacco smoke | never-smoker mean; current-smoker mean | never-smoker ≈ 0.85 |
| cg21566642 | tobacco smoke (secondary) | never-smoker mean | — |
| cg06690548 ★ | SLC7A11 — alcohol intake | means by intake band | light-moderate ≈ 0.86 |
| cg17860381, cg20241083, cg03546163 | NR3C1 / FKBP5 — glucocorticoid signalling | cohort mean + SD | +2.4σ vs cohort |
| cg09649689, cg16232979, cg11761638 | CLOCK / PER2 / ARNTL — circadian | cohort mean + SD | +1.1σ vs cohort |
| 19-site composite | long-term residential PM2.5 | composite score distribution | +0.7σ vs cohort |
| 7-site composite | early-life adversity | composite score distribution | −0.2σ vs cohort |

Composite scores need the site list and the weighting as published, not just a
distribution — record those alongside the reference values.
