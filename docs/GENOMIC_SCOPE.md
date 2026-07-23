# Scoping: genomic data, phenotypes, and how MethylAsk relates to the genomics work

*Grounded in code and data that already exist — the "One Person, Six Genomes" project (Claude Science project `proj_66ffa9d0a3de`, committed to the personal `colbyt/genomics` GitHub repo) — not a hypothetical.*

## The end state we want

One place a citizen uploads a file — a 23andMe export, a VCF, or an Oxford Nanopore run — and gets **all relevant analysis back**: variant/clinical interpretation *and* methylation interpretation, in one report. The worry is repo sprawl and confusion. This document proposes a structure that reaches that end state with the fewest moving parts.

## What already exists (so we don't rebuild it)

The genomics project already built and validated the variant half. Reusable pieces pulled into this project:

| Piece | File | What it does |
|---|---|---|
| 23andMe → VCF | `23andme_to_vcf.py` | Reference-anchored VCF from consumer raw data (handles no-calls, haploid X/Y/MT, indel placeholders) |
| Complete Genomics → VCF | `cg_var_to_vcf.py` | Reconstructs a proper callset from a CG var file |
| Build lift | `lift_b36_posmap.py` | GRCh36/37 → GRCh38 coordinate lift |
| Consensus merge | `build_consensus.py` | Merge multi-source genomes into one provenance-tagged callset |
| Reference-genotype resolver | `refgt_resolver.py` | Confirms hom-reference at sites absent from a callset (the "normal here, not missing" problem) |
| Trait framework | `trait_report.py` | Given a trait table (rsid, effect_allele, trait), reports genotype + multi-source confidence |
| ClinVar panel cache | `clinvar_panel_157genes.json.gz` | 393,806 variants across 157 medically-relevant genes (ACMG SF v3.2 + carrier + pharmacogenes) |
| Unified genome | `cameron_thomson.unified.GRCh38.vcf.gz` | 4.11M variants, per-site provenance (a real end-to-end test input) |

The key realization: **MethylAsk and the genomics work already share the same shape** — provider-style lookups against cached reference databases (ClinVar, gnomAD), coordinate normalization across builds, evidence-tiered findings, multi-source confidence, and a single readable report. The genomics project even validated against a **2018 methylation panel** — it already reached the genome↔methylome seam from the other side.

## The decision that actually matters: variant *interpretation* vs variant *calling*

Two different weights of work hide under "genomic analysis," and they belong in different places:

1. **Interpretation** (light): a VCF or 23andMe file is already variant calls. Reading them is lookups + coordinate joins + evidence tiering — the same profile as MethylAsk. All the pulled code above is interpretation. **This is in scope and cheap.**
2. **Calling** (heavy): raw reads (BAM / modBAM) → a VCF. Alignment, a variant caller, filtering, QC — a real bioinformatics pipeline with its own tool stack (minimap2/dorado, Clair3/DeepVariant, bcftools), real CPU-hours, and a higher clinical-liability surface. **This is the only genuinely different animal.**

Most DIY-citizen inputs are already-called (23andMe, an uploaded VCF), so interpretation covers the common case. Calling matters specifically for the ONT case below.

## The ONT case is why a shared core wins

Oxford Nanopore produces a **modBAM**: the *same reads* carry both the aligned sequence (→ variants) and per-base methylation (→ methylome), in one file from one run. They are physically inseparable at the source. MethylAsk's design already lists modBAM as a phase-2 input.

This forces the architecture: **the modBAM reader must live in a shared layer**, parse once, and emit two streams — methylation to the methylation engine, reads/variants to the variant engine. If methylation and genomics were fully separate silos, that file would be parsed twice by two codebases with two coordinate conventions. If they were one monolith, MethylAsk's clean lightweight footprint would absorb the whole variant-calling stack.

## Recommendation — three repos, one shared core, minimal sprawl

Sprawl comes from *duplicated or unclear-boundary* repos, not from repo *count*. Three repos with crisp, non-overlapping responsibilities are clearer than two that both grow toward each other. Proposed:

```
  goodancestor/bio-core        (shared library — the anti-sprawl keystone)
     ├─ provider interface, cache, error-tolerant status
     ├─ evidence tiering (robust/moderate/speculative/unknown)
     ├─ coordinate normalization (build lift, probe/rsID → locus)
     ├─ single-disclaimer model + report renderer (HTML/PDF)
     ├─ ClinVar / gnomAD / dbSNP provider (shared by both)
     └─ modBAM reader → (methylation stream, read stream)

  goodancestor/methylask       depends on bio-core
     └─ methylation ingestion, normalization, epigenetic clocks,
        EWAS/GDC/methylation-DB interpretation

  goodancestor/genomics        depends on bio-core   (generalized from colbyt/genomics)
     └─ VCF / 23andMe / CG ingestion, multi-source consensus,
        trait framework, ClinVar clinical screen, (optional) variant calling
```

Then the "combined" experience is **not a fourth codebase** — it is a thin application that imports both engines:

```
  goodancestor/report          depends on methylask + genomics + bio-core
     └─ upload one file → detect type → route to the right engine(s) →
        one merged, evidence-tiered report
```

`report` is the "simple upload, all relevant analysis" front door. It stays thin because all the real work lives in the two engines and the shared core. Detection logic is small: `.txt`/23andMe → genomics; IDAT/β-matrix → methylask; VCF → genomics; **modBAM → bio-core splits it → both**.

### Why this beats the alternatives

- **vs. one merged monolith:** keeps MethylAsk's light footprint and separate release cadence; the heavy variant-calling stack only loads in `genomics`. Methylation users don't pay for a variant-calling dependency tree.
- **vs. two disconnected silos:** no duplication of the provider interface, tiering, disclaimer, coordinate normalization, or the ClinVar cache — and the ONT file is parsed once. This is the actual sprawl-and-confusion risk, and the shared core removes it.
- **vs. keeping genomics personal (`colbyt/genomics`):** generalizing it to `goodancestor/genomics` de-personalizes the code (your data becomes a test fixture, not the subject), which is the prerequisite for the eventual public release.

### Migration path (low-risk, incremental)

1. **Lift, don't rewrite.** The pulled scripts already work. First move is to relocate the truly-shared pieces (ClinVar provider + cache, coordinate normalization, trait framework) into `bio-core`, leaving thin wrappers behind.
2. **Generalize `genomics` in place.** Replace hard-coded personal paths/sample names with parameters; keep `cameron_thomson.unified.GRCh38.vcf.gz` as an integration-test fixture.
3. **MethylAsk adopts `bio-core`** for its ClinVar provider and normalization instead of its own copies (they were built to match anyway).
4. **`report` last** — only once both engines expose a stable "analyze(file) → findings" call.

All repos private until reviewed; public release is a later, gated step and a non-goal for now.

## Where plants fit — split along knowledge, not organism

The seagrass (plant) epigenomics work raises a fair question: is there credible overlap we'd waste by ignoring plants, or do they need separate plant-vs-animal sections? The resolution is to split along the right axis.

**Two kinds of code live in all of this:**
- **Mechanism** — how you read and move data: parse a modBAM, pile up a bedMethyl, normalize a VCF, intersect features, lift coordinates, render a report. **Organism-agnostic** — a human modBAM and an eelgrass modBAM are the same format read by the same code.
- **Knowledge** — what a marker means: ClinVar pathogenicity, Horvath's clock CpGs, the RdDM pathway's CHH targets. **Deeply organism-specific**, transfers not at all.

`bio-core` is **pure mechanism, and therefore has no species.** That is where the plant/human overlap actually lives. Organism *knowledge* stays in its own repo:

```
  goodancestor/bio-core     mechanism, organism-agnostic, CONTEXT-AWARE
  goodancestor/methylask    human methylation knowledge (EWAS, clocks, Illumina arrays)
  goodancestor/genomics     human variant knowledge (ClinVar, PRS, pharmacogenomics)
  goodancestor/seagrass     plant epigenomics knowledge (RdDM, TE methylation, breeding)
        ↑ all depend on bio-core; none duplicate the mechanism
```

No plant-vs-animal *sections* inside a package (that is sprawl inside a repo). Plant knowledge lives in `seagrass`, but `seagrass` **depends on bio-core** for the plumbing instead of carrying its own copy — the difference between "ignore plant people" (wasteful) and "share the plumbing, not the meaning" (right). Arrays and clinical DBs stay human because they have no plant analog — they correctly live in methylask/genomics, never in bio-core.

### The one design choice that makes bio-core plant-ready — at zero cost

The load-bearing plant/human divergence is **sequence context**. Human methylation is ~entirely CpG (arrays and clocks assume it). Plants use CpG + CHG + CHH, and CHH (RdDM-deposited) is the *central* signal in the seagrass work, not incidental. So bio-core's methylation data model must represent *"a cytosine, its context (CG/CHG/CHH), and its level"* — not *"a CpG and its beta."* This is strictly more correct even for humans (non-CpG methylation in brain/stem cells is real), and both bedMethyl and modBAM already carry the context code. Building bio-core context-aware is not a plant concession; it is the correct data model, and it keeps the plant door open for free.

## Tooling audit — what the seagrass project actually exercised

Audited the seagrass project's methods appendix and code artifacts for liftable mechanism. Finding: **no custom seagrass analysis code to pull** (the methylation landscape ran as a streaming awk pass; everything else is standard tools invoked directly), but a **validated specification** of exactly which mechanical primitives bio-core needs — a second independent project confirming they are organism-agnostic.

| Mechanism seagrass used | Detail | To bio-core? |
|---|---|---|
| bedMethyl parsing | modkit 18-col, context in col4 (`m,CG,0`); weighted methylation = ΣN_mod / Σ(N_mod+N_canonical) at cov≥5× | **Yes — verbatim.** Most reusable piece; MethylAsk needs identical logic |
| Sequence-context partition | per-context CG/CHG/CHH weighted methylation | **Yes** — empirically justifies the context-aware model above |
| VCF merge + biallelic filter | `bcftools merge -0`, `view -m2 -M2 -v snps` | **Yes** — same bcftools wrappers genomics uses |
| Feature intersection | `bedtools intersect -u -sorted` over genes/exons/introns/promoters/repeats | **Yes** — MethylAsk's probe→feature annotation is the same op |
| Coordinate concordance QC | verify positions fall within contig bounds | **Yes** — belongs in bio-core normalization/QC |
| Genetic distance / Mantel / PCA | scikit-allel IBS, Ward linkage, Mantel | Partly — population genetics; leans genomics, not core |
| DIAMOND ortholog search | `diamond blastp --very-sensitive` | No — comparative genomics, stays in seagrass |
| De novo TE annotation | RepeatModeler / RepeatMasker / EDTA | No — plant-specific, stays in seagrass |

Two concrete design consequences:
1. **bio-core's bedMethyl reader implements the modkit 18-column contract and the cov≥5×, ΣN_mod/Σ(N_mod+N_canonical) weighted-methylation formula verbatim** — both projects need exactly this.
2. **The context-aware data model is empirically justified, not speculative** — seagrass's headline result (CHH/RdDM) cannot be represented in a CpG-only model.

*Caveat:* only 2 frames and 7 code artifacts were reachable in the shared project from here; deeper seagrass sessions and any `GoodAncestor/seagrass` repo were not accessible with the current token, so this audit rests on the methods appendix rather than raw source. A token scoped to the seagrass repo would let us confirm whether any wrapper code there is worth lifting directly.

## Two open decisions (flagged, not blocking)

1. **Is variant *calling* in scope, or interpretation only?** If citizens arrive with VCF/23andMe files, `genomics` stays light and calling can wait. If raw ONT read → variant calling is a launch feature, that pulls in the heavy stack and is the main scoping cost.
2. **PRS and pharmacogenomics — in or out?** Highest consumer value, highest clinical-claim sensitivity. The existing work already surfaces pharmacogenomic findings (CYP4F2, DPYD, NAT2…) at a research-information framing, which fits the "no definitive health claims" stance — but PRS in particular needs careful evidence-tier labeling before it faces users.
