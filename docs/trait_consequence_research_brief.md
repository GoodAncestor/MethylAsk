# Research brief — trait consequence copy for methylation reports

Hand this to a research agent. Self-contained; paste it whole.

---

## Context

A methylation report currently tells a reader that a site on their DNA is
"associated with lower methylation" for a given trait, and shows their own
reading. What it cannot yet tell them is **what that means and why they should
care**. This brief fills that gap.

The traits come from the EWAS Catalog (8,023,135 published associations across
6,515 free-text trait strings). We have already computed the frequencies, so the
list below is not a guess — it is what actually appears in reports, ranked.

## The distinction that matters most

**An EWAS association is a population-level statistical finding. It is not a
diagnosis, not a risk score, and not a prediction about the individual reading
the report.** Copy that blurs this is the single failure mode that matters here.

- ✅ "Studies comparing people with and without rheumatoid arthritis have found
  differences in methylation at sites like this one."
- ❌ "This marker indicates rheumatoid arthritis risk."
- ❌ "Your methylation pattern suggests you may develop this condition."

Most EWAS findings are **cross-sectional**: methylation and trait were measured
at the same time, so the association cannot say which came first. For several
traits the methylation change is widely believed to be a *consequence* of the
trait, not a cause of it — smoking is the clearest example. Where the literature
has a view on direction of causation, say so; where it does not, say that.

Never write copy that implies the reader has, or will get, the condition.

## What we need, per trait

Four short prose fields, each 1-3 sentences, plainly written for a general adult
reader with no biology background:

1. **`what_it_is`** — what the trait is. Assume nothing. "HDL cholesterol" needs
   a sentence saying what HDL is.
2. **`what_an_association_means`** — what it means that methylation at some site
   is associated with this trait. This is the "why should I care" field, and it
   is where the study design belongs: is this typically case-control, a
   continuous measure in a cohort, prospective or cross-sectional?
3. **`what_it_is_not`** — the limits. What a reader must not conclude. Be
   specific to this trait rather than reusing a generic disclaimer.
4. **`typical_evidence`** — one sentence on how well replicated this trait's
   methylation associations are overall. Smoking is among the most replicated
   findings in the field; many others rest on a single cohort. Say which.

Plus `pmids` — the sources these statements rest on. Reviews and consortium
papers are ideal here; this is background, not a numeric claim.

## Source quality

Peer-reviewed, human, methylation-specific where possible. Prefer review
articles and consortium/meta-analysis papers over single cohorts for background
statements. If a trait's methylation literature is thin, that is a finding —
record it in `typical_evidence` rather than padding.

## Rules that are not negotiable

1. **No medical advice, no risk language, no actionability.** Do not tell a
   reader what to do, and do not say whether a higher or lower reading is better.
2. **No valence.** Do not label a trait or a direction as good or bad. State
   what was measured. (High HDL is conventionally called "good cholesterol" —
   even that framing is a judgement; describe what HDL does instead.)
3. **Never imply the reader has the condition.** See above; this is the one that
   would do real harm.
4. **Report "thin literature" as a real result.** A trait with little
   methylation-specific work is expected and handled.
5. **Every substantive claim carries a PMID.**

## Traits wanted

Ranked by how many published associations carry them, so this is exactly the
order of impact. Percentages are of all 8,023,135 associations.

| Trait | Associations | % |
|---|---|---|
| age | 690,552 | 8.61% |
| BMI / body mass index | 158,721 | 1.98% |
| HDL cholesterol | 134,478 | 1.68% |
| body fat | 105,353 | 1.31% |
| C-reactive protein (CRP) levels | 95,492 | 1.19% |
| sex | 81,547 | 1.02% |
| waist-hip ratio | 78,737 | 0.98% |
| HIV infection | 67,506 | 0.84% |
| gestational age | 65,962 | 0.82% |
| COPD (incident) | 55,873 | 0.70% |
| rheumatoid arthritis | 52,367 | 0.65% |
| smoking / tobacco smoking | 40,124 | 0.50% |
| clear cell renal carcinoma | 35,550 | 0.44% |
| type 2 diabetes | 32,922 | 0.41% |
| eosinophilia | 24,114 | 0.30% |
| total cholesterol | 19,378 | 0.24% |
| alcohol consumption | 14,304 | 0.18% |
| atopy | 13,556 | 0.17% |
| primary Sjögren's syndrome | 12,238 | 0.15% |
| pancreatic ductal adenocarcinoma | 11,634 | 0.15% |
| schizophrenia | 10,173 | 0.13% |
| chronic pain | 8,651 | 0.11% |

Casing and phrasing variants ("BMI"/"Body mass index", "sex"/"Sex") are the same
concept — write one entry per concept.

**Two categories deliberately excluded**, for information rather than action:
study covariates such as "Tissue" (328,363 associations — a study design
variable, not a property of a person), and protein-level measurements, which are
handled separately below.

## The protein group — one template, not 130 entries

**64% of all associations are protein abundance measurements**, appearing either
as bare UniProt accessions (P02748, O75882) or as "PAPPA protein levels
(SeqId = 4148-49)". Do NOT write per-protein copy; the individual protein names
and functions are being resolved separately from UniProt.

What is needed instead is **one carefully written explanation of the whole
class**: what it means that methylation at a site is associated with the blood
level of some protein, why so much of this literature exists (large proteomic
panels run against methylation arrays), and what a reader can and cannot take
from it. This one piece of copy carries two thirds of the report, so it deserves
disproportionate care.

## Output format

**1. A JSON object**, loaded directly:

```json
{
  "smoking": {
    "label": "Tobacco smoking",
    "what_it_is": "...",
    "what_an_association_means": "...",
    "what_it_is_not": "...",
    "typical_evidence": "...",
    "pmids": ["23691101"]
  },
  "_protein_level": {
    "label": "Blood level of a protein",
    "what_it_is": "...",
    "what_an_association_means": "...",
    "what_it_is_not": "...",
    "typical_evidence": "...",
    "pmids": ["..."]
  }
}
```

Keys are lowercase canonical concept names. `_protein_level` is the class-wide
entry described above. Omit any trait you could not source — do not emit a
placeholder.

**2. A provenance table** — one row per trait: concept, the PMIDs used, first
author + year for each, and what kind of source it was (review, meta-analysis,
single cohort). Plus an explicit list of traits where the methylation literature
was too thin to support the four fields, and what you searched.

## Sanity checks before returning

- Read each `what_it_is_not` as if you were a worried person who just saw their
  own number next to it. If any sentence could be read as "you might have this",
  rewrite it.
- No field should contain the words "risk of", "indicates", "suggests you", or
  any recommendation.
- Confirm every PMID resolves to the paper you cite.
