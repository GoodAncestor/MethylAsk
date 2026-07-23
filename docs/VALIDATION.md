# MethylAsk — Data Source Validation Results

*Live probe of every candidate data source, run 2026-07-23. Each endpoint was contacted directly and its response inspected. Status values feed the "error-tolerant provider" behavior in the design: a source that errors is kept and labeled, not dropped.*

## Working — validated with real methylation/clinical data

| Source | Endpoint | Result | Access notes |
|---|---|---|---|
| **EWAS Catalog** | `www.ewascatalog.org/api/?cpg=<probe>` | `cg00000029` → 12 associations (trait, gene, chrpos, beta, p, n, tissue, array, PMID) | JSON; rows are **positional arrays** paired with a top-level `fields` list — zip them. Required an allowlist grant. |
| **GDC / TCGA** | `api.gdc.cancer.gov/files` | 20,397 "Methylation Beta Value" files; SeSAMe level-3 β, EPIC/450K | Rich filter+field query language. Required an allowlist grant. |
| **ClinVar** | NCBI E-utilities | `BRCA1[gene]` → 16,041 variant records | Allowlisted by default. |
| **dbSNP** | NCBI E-utilities | esummary returns rsID records | Allowlisted by default. |
| **GWAS Catalog** | `www.ebi.ac.uk/gwas/rest` | `rs7329174` → functionalClass `intron_variant` | Allowlisted by default. |
| **Ensembl** | `rest.ensembl.org/overlap` | region query → 7 gene features | Allowlisted; one call timed out then succeeded (rate/latency). |
| **UCSC** | `api.genome.ucsc.edu` | genome list OK | Allowlisted. |
| **openFDA** | `api.fda.gov/drug/label` | label records OK | Allowlisted; drug-label layer for pharmacogenomics. |
| **ENCODE** | `www.encodeproject.org/search` | WGBS experiment search OK | Allowlisted; regulatory-context layer. |
| **Zhou-lab InfiniumAnnotation** | GitHub API + release assets | repo reachable | Allowlisted; bulk manifest/annotation download. |
| **Illumina** | `support.illumina.com` EPICv2 product files | page reachable | Required an allowlist grant; manifest download. |

## Errored — kept in design, labeled, retry later

| Source | Endpoint | Error | Nature | Handling |
|---|---|---|---|---|
| **EWAS Atlas / EWAS Open Platform** | `ngdc.cncb.ac.cn/ewas/api` | `502 Bad Gateway` (persistent across paths) | **Server-side outage**, not an access block (allowlist grant succeeded) | Provider registered; `status()` reports `unavailable: upstream 502`. EWAS Catalog covers overlapping data meanwhile. Retry on schedule. |
| **OpenGWAS (MRC-IEU)** | `gwas-api.mrcieu.ac.uk` | `SSL: CERTIFICATE_VERIFY_FAILED — certificate expired` | **Server-side cert expiry** on their host | Provider registered; `status()` reports `unavailable: expired TLS cert`. Do **not** disable cert verification as a workaround. Retry on schedule; flat-file download is the fallback ingest path. |

## Download sizes — where each source lives (GitHub vs NAS)

Measured live from HTTP headers on 2026-07-23. **Storage capacity is not the constraint** — everything, including the 294 GB GDC corpus, fits on a normal NVMe drive on the hosting server. The only hard rule is GitHub's 100 MB per-file limit for what gets committed to the repo; everything else lives on the server's local disk. Rule of thumb: **≤~100 MB and static → commit to the GitHub repo; larger or frequently rebuilt → server disk**.

| Source (flat file) | Compressed size | Home | Why |
|---|---|---|---|
| Zhou HM450 hg38 manifest | 22 MB | **GitHub** | Static reference, small |
| Zhou EPIC (850K) hg38 manifest | 37 MB | **GitHub** | Static reference |
| Zhou EPICv2 hg38 manifest | 39 MB | **GitHub** | Static reference |
| EWAS Catalog studies table | 0.27 MB | **GitHub** | Tiny metadata |
| Epigenetic clock coefficients | <5 MB (est.) | **GitHub** | Small, static, ships with tool |
| EWAS Catalog full results dump | 174 MB | **Server disk** | Over the git 100 MB limit; rebuilt on release |
| ClinVar GRCh38 VCF.gz | 193 MB | **Server disk** | Over limit; updated weekly by NCBI |
| GDC/TCGA methylation β corpus | **294 GB** (20,397 files, avg 14.4 MB) | **Server disk (full mirror)** | Too large for git; fits on NVMe and downloads in <1 h on 10 Gbit, so mirror the whole corpus rather than cache on demand |
| OpenGWAS trait dumps | GB-scale per dataset | **Server disk** | Bulk flat files; the chosen path for OpenGWAS |

Practical guidance: GitHub rejects single files >100 MB (Git LFS needed beyond that). Everything tagged **GitHub** fits comfortably as normal committed files. Everything tagged **Server disk** either exceeds the 100 MB single-file limit or changes often enough that a synced cache on local disk is cleaner than versioning it in git. The 80 TB NAS is available for backup/mirroring but is not required for capacity.

## OpenGWAS decision (prototype)

The expired-certificate question is moot for the build: OpenGWAS is registered as a **download-backed synced cache** on the NAS, not a live API call. This removes both the flapping endpoint and the TLS question from the runtime path, fits the "keep up on a schedule" model, and is appropriate because the data is read-only public trait associations. If a live query is ever wanted during prototyping, TLS verification can be disabled for that one provider only (read-only public data, no PHI on the wire) — but it is not needed for the cache-backed path.

## Takeaways for the build

- The interpretation layer is viable today: EWAS Catalog + GDC + ClinVar + GWAS Catalog + ENCODE all return usable per-marker data live.
- Two sources are temporarily down on their own infrastructure. Under the error-tolerant provider model they stay in the registry with a recorded status and a scheduled retry, and their reports carry a "source unavailable at generation time" note rather than silently omitting them.
- Every allowlist grant needed for a private-cloud deployment is now known and documented (see design §6).
