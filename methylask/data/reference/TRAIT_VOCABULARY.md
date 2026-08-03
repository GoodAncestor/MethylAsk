# Trait vocabulary — `trait_canonical.json` and `protein_traits.json`

Derived reference data describing the EWAS Catalog's free-text trait vocabulary.
Both were built from the catalog's public bulk downloads on 2026-08-03, over
**8,023,135 associations across 6,515 distinct trait strings**.

Neither file is a source of scientific claims. They classify and name what the
catalog already contains, so a report can decide what to show and what to call it.

## `trait_canonical.json`

The top 400 non-accession traits by association count, clustered into 381
canonical concepts. Per entry: `canonical_label`, `variants` (every raw string
that maps to it), `total_associations`, `class`.

Verified: variant counts sum to each concept total, every raw string maps exactly
once, and the file totals 2,549,800 — matching the true top-400 sum computed
independently from the catalog.

### Classes, and why `covariate` matters

| class | concepts | associations | % of catalog |
|---|---:|---:|---:|
| `health_trait` | 164 | 1,947,864 | 24.28% |
| `covariate` | 4 | 337,194 | 4.20% |
| `other` | 10 | 189,045 | 2.36% |
| `protein_level` | 203 | 75,697 | 0.94% |

**`covariate` entries are study-design variables, not properties of a person**,
and should not be presented as findings about the reader:

| concept | associations |
|---|---:|
| Tissue | 328,363 |
| Age × Sex interaction term | 6,014 |
| Maternal education at time of pregnancy | 2,390 |
| Ethnicity | 427 |

`Tissue` alone is the second-largest trait in the entire catalog. Rendering it as
"Tissue — associated with lower methylation at this site" is meaningless to a reader.

`other` holds concepts that could not be confidently placed — `Sex`, `Infant sex`,
`Fetal vs adult liver`, `Age 4 vs age 0`, and four Metabolon metabolite readouts.
For display purposes these behave like covariates even though the class differs.

Two concepts are flagged as genuinely ambiguous and deliberately left in
`health_trait`: **Age** (717,205 — the largest concept in the dataset) and
**Gestational age**. Both are legitimate subjects of epigenetic-clock research
*and* near-universal adjustment covariates in EWAS models. The raw strings give
no way to tell which role dominates. Do not resolve this silently.

Merge policy was conservative. Case, spacing and typo variants were merged;
`Incident X` / `Prevalent X (Self-report)` / bare `X` were kept separate across
~13 diseases, as were maternal/paternal variants of the same measure.

## `protein_traits.json`

Protein-abundance traits — **5,514,449 associations, 68.73% of the entire
catalog** — as 4,131 distinct proteins. They appear either as bare UniProt
accessions (`P02748`) or written as `PAPPA protein levels (SeqId = 4148-49)`.

Per entry: `raw_trait_strings`, `accession`, `gene_symbol`, `protein_name`,
`function_sentence`, `total_associations`.

| | count | associations |
|---|---:|---:|
| with accession | 4,118 | — |
| with a UniProt function sentence | 3,877 | 4,850,035 (88% of the group) |
| no accession | 13 | 668 |
| no function annotation | 254 | 664,414 |

Resolved by joining against the reviewed human proteome
(`organism_id:9606 AND reviewed:true`, 20,431 entries) rather than per-protein
API calls. `function_sentence` is the **first sentence** of UniProt's own
`Function [CC]` annotation, with `{ECO:...}` evidence codes and inline
`(PubMed:...)` references stripped. It is quoted description, not interpretation —
it never says whether a level is good or bad.

The four largest unresolved proteins (`LRG1`, `A1BG`, `LUM`, `A0A0G2JMB2`) have
accessions but no UniProt function annotation at all. That is a genuine gap in
the source, not a fetch failure — do not fill it from elsewhere without saying so.

### Attribution

`protein_traits.json` contains text from **UniProt, licensed CC BY 4.0** —
attribution is a licence condition. Registered in bio-core as the `uniprot`
source so it appears in the report's Data sources panel.

## Regenerating

Both derive from:

    https://www.ewascatalog.org/static/docs/ewascatalog-results.txt.gz   (~174 MB)
    https://www.ewascatalog.org/static/docs/ewascatalog-studies.txt.gz   (~272 KB)

joined on `study_id` to attribute each association to its study's trait, then
counted. Protein resolution additionally uses the UniProt REST `stream` endpoint.

Note for anyone re-deriving the accession split: a UniProt accession pattern of
`^[A-Za-z][A-Za-z0-9]{5,9}$` **is wrong** — it matches ordinary words like
`Tissue`, `Smoking`, `Melanoma` and `Asthma`, silently dropping real traits.
Use `^[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z0-9]{3}[0-9])?$`.
