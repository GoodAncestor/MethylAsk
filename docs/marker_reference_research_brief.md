# Research brief — populate `marker_reference.json`

Hand this to a research agent. It is written to be self-contained; paste it whole.

---

## What we need

Absolute **population mean beta values** for a fixed list of CpG probes, each in a
named reference group, quoted from published research with the PMID recorded.

These feed a personal methylation report that tells someone "your value at
cg05575921 is 0.71; never-smokers average 0.85 [source]". The number is a claim
about a population and it must belong to a paper, not to us.

## The one distinction that matters most

**We need absolute levels, not effect sizes.** Most EWAS publications report the
*difference* a trait makes — a Δβ, a regression coefficient, a log-odds, a
percentage-point change. Those are NOT what this table holds and cannot be
converted into it without the baseline they were measured against.

- ✅ "Mean beta at cg05575921 was 0.85 (SD 0.03) in never-smokers and 0.68 in
  current smokers" — usable, both groups.
- ❌ "Smoking was associated with a −0.17 change in beta at cg05575921
  (p = 3×10⁻⁴⁵)" — an effect size. Not usable on its own.
- ❌ "cg05575921 was hypomethylated in smokers" — directional only. Not usable.

If a paper gives an effect size *and* a stated reference-group mean, take the
mean and note the effect size separately. If it only gives the effect size, that
paper does not answer this question — keep looking or report nothing found.

The values usually live in a descriptive-statistics table, a supplementary
table, or the axis/annotation of a distribution figure — rarely in the abstract.

## Source quality

Required: peer-reviewed; human; reports the value for a **named group**
(never-smoker, current smoker, cohort mean, intake band); states the **tissue**
and the **array platform**; states the **n** for that group.

Prefer, in order: meta-analyses and consortium studies (large n) > single large
cohorts > small cohorts. Record the best available, not the first found.

Tissue is not interchangeable — a whole-blood value does not transfer to buccal
or saliva. Record the tissue exactly as the paper states it. Array platform
(27k / 450k / EPIC / EPICv2) should be recorded even though beta values are
broadly comparable across them.

## Rules that are not negotiable

1. **Never estimate, interpolate, or infer a number.** If the paper does not
   print it, it does not exist for our purposes.
2. **`sd` is optional — omit it unless the paper published a standard deviation
   for that group.** Do not derive an SD from a confidence interval, standard
   error, IQR, or figure. A downstream sigma computed from a guessed spread is a
   fabricated statistic on a health report.
3. **Report "not found" as a real result.** A marker with no published reference
   value is expected and handled — the report says so. An invented value is far
   worse than a gap.
4. **Every entry needs a PMID and a locator** — the table, figure, or supplement
   number where the value appears, so it can be spot-checked in under a minute.

## Markers wanted

| Probe(s) | Signature | Reference groups wanted |
|---|---|---|
| cg05575921 | AHRR — tobacco smoke | never-smoker, former smoker, current smoker |
| cg21566642 | tobacco smoke (secondary marker) | never-smoker, current smoker |
| cg06690548 | SLC7A11 — alcohol intake | means by intake band (none / light-moderate / heavy) |
| cg17860381, cg20241083, cg03546163 | NR3C1 / FKBP5 — glucocorticoid signalling | general-population cohort mean + SD |
| cg09649689, cg16232979, cg11761638 | CLOCK / PER2 / ARNTL — circadian | cohort mean + SD; day-shift vs night-shift if published |
| 19-site composite | long-term residential PM2.5 | see composites below |
| 7-site composite | early-life adversity | see composites below |

### Composites

For the two composite scores, the distribution alone is not enough. Return
**the published probe list and the weighting or scoring method**, plus the score
distribution (mean + SD) in the population it was derived in. Without the site
list and weights the score cannot be reproduced, and a distribution for a score
we compute differently is meaningless. If several competing composites exist,
report each separately rather than merging them.

## Output format

Return two things.

**1. A JSON object** in exactly this shape (this is the file schema — it is
loaded directly):

```json
{
  "cg05575921": {
    "gene": "AHRR",
    "references": [
      {
        "group":  "never-smoker",
        "beta":   0.85,
        "sd":     0.03,
        "tissue": "whole blood",
        "array":  "450k",
        "n":      1793,
        "pmid":   "23691101"
      }
    ]
  }
}
```

Omit `sd` entirely when unpublished. Multiple groups per marker go in the same
`references` list. Omit any marker for which nothing was found — do not emit a
placeholder entry.

**2. A provenance table**, one row per value returned: probe, group, value,
PMID, first author + year, and the exact locator (e.g. "Supplementary Table 3",
"Table 1, row 4", "Figure 2B axis"). Plus a short list of markers searched with
nothing usable found, and what was searched, so the gap is documented rather
than silent.

## Sanity checks before returning

- Beta values are proportions and must fall in 0–1. A number like 85 or −0.17 is
  a percentage or an effect size, not a beta — recheck it.
- A never-smoker AHRR value near 0.85 and a current-smoker value clearly lower is
  the expected shape; if the direction is inverted, confirm the paper's beta
  convention before returning it.
- Cross-check any headline value against a second publication where one exists,
  and note agreement or disagreement.
