"""Methylation findings are explained as group biology, with the reading in context."""
from biocore.providers.base import Finding, Tier, Category
from methylask.meaning import interpret


def test_ewas_uses_the_trait_copy_and_the_reading():
    f = Finding(marker="cg00574958", source="ewas_catalog", description="x", tier=Tier.ROBUST,
                categories=[Category.CLINICAL],
                detail={"trait": "Body mass index", "copy_key": "bmi", "beta": 0.02, "p": 1e-12,
                        "n": 5387, "tissue": "whole blood", "your reading": 0.41, "topic": "metabolic"},
                pmids=["28213390"])
    assert interpret([f]) == 1
    ip = f.interpretation
    assert "BMI" in ip.found and "0.41" in ip.found
    assert "not a measurement of the reader" in ip.can_mean.lower()
    assert "5,387" in ip.how_sure and "whole blood" in ip.how_sure
    assert ip.next_step == "" and ip.copy_version == "trait_copy"
    assert [c.kind for c in f.evidence_chain][:2] == ["variant", "trait"]


def test_ewas_without_copy_falls_back_to_the_generic_sentence():
    f = Finding(marker="cg1", source="ewas_catalog", description="x", tier=Tier.MODERATE,
                categories=[Category.TRAIT], detail={"trait": "Height", "p": "bad", "n": None})
    interpret([f])
    assert "not a measurement of your trait" in f.interpretation.can_mean
    assert f.interpretation.how_sure == "" and f.interpretation.copy_version == "inline"


def test_gdc_is_reference_biology_not_a_screen():
    f = Finding(marker="cg1", source="gdc", description="x", tier=Tier.ROBUST,
                categories=[Category.CLINICAL],
                detail={"project": "TCGA-BRCA", "delta_beta": 0.31, "n_tumor": 780, "n_normal": 96,
                        "your reading": 0.22, "topic": "cancer"},
                link="https://portal.gdc.cancer.gov/projects/TCGA-BRCA")
    interpret([f]); ip = f.interpretation
    assert "tumour tissue" in ip.found.lower() and "TCGA-BRCA" in ip.found and "+0.31" in ip.found
    assert "not a cancer test" in ip.can_mean.lower()
    assert "780" in ip.how_sure and "96" in ip.how_sure


def test_other_sources_are_left_alone():
    f = Finding(marker="cg1", source="clocks", description="x", tier=Tier.UNKNOWN, categories=[Category.AGING])
    assert interpret([f]) == 0 and f.interpretation is None
