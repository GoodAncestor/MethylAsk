# Epigenetic clock coefficients — provenance

Published coefficients extracted verbatim from the `dnaMethyAge` R package
(github.com/yiluyucheng/dnaMethyAge, `data/*.rda`), which curates them from the
original papers. Each CSV is `Probe, Coefficient`; where present, an `Intercept`
row carries the model intercept.

| File | Clock | CpGs | Intercept | Age transform | Origin |
|---|---|---|---|---|---|
| Horvath2013_PanTissue.csv | Horvath pan-tissue | 353 | yes | anti-log-linear (adult age 20) | Horvath 2013, Genome Biol |
| Hannum2013_Blood.csv | Hannum blood | 71 | no | linear | Hannum et al. 2013, Mol Cell |
| Levine2018_PhenoAge.csv | PhenoAge | 513 | yes | linear | Levine et al. 2018, Aging |
| Horvath2018_SkinBlood.csv | Skin & Blood | 391 | yes | anti-log-linear (adult age 20) | Horvath et al. 2018, Aging |

Verify counts: Horvath = 353 CpG (as published), Hannum = 71, PhenoAge = 513,
Skin&Blood = 391.

Note (docs/DESIGN.md §3c): some clock CpGs are absent from EPIC/EPICv2 arrays.
The engine reports coverage (how many of the clock's CpGs were found in the
sample) so a low-coverage estimate can be flagged rather than silently biased.
