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
    assert ip.found == "Body mass index — methylation here rises with it."
    assert f.detail["short_label"] == "Body mass index"
    assert "Not a measurement of you" in ip.can_mean
    assert ip.how_sure == "" and ip.next_step == "" and ip.copy_version == "trait_copy"
    assert [c.kind for c in f.evidence_chain][:2] == ["variant", "trait"]


def test_ewas_without_copy_falls_back_to_the_generic_sentence():
    f = Finding(marker="cg1", source="ewas_catalog", description="x", tier=Tier.MODERATE,
                categories=[Category.TRAIT], detail={"trait": "Height", "p": "bad", "n": None})
    interpret([f])
    assert f.interpretation.found == "Height — methylation here differs with it."
    assert f.interpretation.copy_version == "inline"


def test_gdc_is_reference_biology_not_a_screen():
    f = Finding(marker="cg1", source="gdc", description="x", tier=Tier.ROBUST,
                categories=[Category.CLINICAL],
                detail={"project": "TCGA-BRCA", "delta_beta": 0.31, "n_tumor": 780, "n_normal": 96,
                        "your reading": 0.22, "topic": "cancer"},
                link="https://portal.gdc.cancer.gov/projects/TCGA-BRCA")
    interpret([f]); ip = f.interpretation
    assert ip.found == "Differs in TCGA-BRCA tumour tissue (+0.31)."
    assert "Not a cancer test" in ip.can_mean and f.detail["short_label"] == "BRCA"


def test_other_sources_are_left_alone():
    f = Finding(marker="cg1", source="clocks", description="x", tier=Tier.UNKNOWN, categories=[Category.AGING])
    assert interpret([f]) == 0 and f.interpretation is None


def test_protein_rows_name_the_protein_not_a_protein():
    f = Finding(marker="cg1", source="ewas_catalog", description="x", tier=Tier.ROBUST,
                categories=[Category.TRAIT],
                detail={"trait": "Blood level of a protein", "protein": "P02748",
                        "subject": "blood level of protein Alpha-2-macroglobulin", "p": 1e-9, "n": 100})
    interpret([f])
    assert f.interpretation.found.startswith("Alpha-2-macroglobulin (a blood protein) — ")


def test_aggregated_row_says_the_replication():
    f = Finding(
        marker="cg1",
        source="ewas_catalog",
        description="x",
        tier=Tier.ROBUST,
        categories=[Category.CLINICAL],
        detail={
            "trait": "Body mass index",
            "copy_key": "bmi",
            "p": 1e-12,
            "n": 9587,
            "n_studies": 3,
            "n_participants": 9587,
            "direction": "consistent",
            "tissues": ["cord blood", "whole blood"],
            "tissue_supported": True,
            "your reading": 0.41,
        },
    )
    interpret([f])
    ip = f.interpretation
    assert ip.found == "Body mass index — methylation here differs with it."
    assert ip.how_sure == ""            # counts and tissues render as chips, not prose
    assert f.detail["n_studies"] == 3 and f.detail["n_participants"] == 9587
