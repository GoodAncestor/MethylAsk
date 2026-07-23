# Committed reference files (<100 MB each)

Small static files that ship in the repo: array manifests (Zhou-lab HM450/EPIC/
EPICv2), epigenetic clock coefficients, EWAS Catalog studies metadata, and the
curated marker -> plain-language table.

Large corpora (GDC 294 GB, ClinVar VCF, EWAS Catalog full dump) are NOT here —
`methylask refresh` mirrors them to data/mirror/ on server disk. See
docs/DESIGN.md §3.4.
