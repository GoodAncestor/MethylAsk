"""Trait class lookup — which findings describe the person vs the study design."""
from methylask.traits import trait_class
from methylask.providers.ewas_catalog import EwasCatalogProvider


def _row(trait, beta=-0.003):
    return {"trait": trait, "gene": "X", "beta": beta, "se": 0.01, "p": 1e-9,
            "n": 100, "tissue": "Whole blood", "methylation_array": "450k",
            "chrpos": "1:1", "pmid": "1"}


def test_tissue_is_a_covariate():
    # 328,363 associations — the catalog's second-largest trait, and a property
    # of the sample rather than of the person
    assert trait_class("Tissue") == "covariate"


def test_a_real_health_trait_is_not_a_covariate():
    assert trait_class("BMI") == "health_trait"


def test_lookup_matches_the_raw_variant_spelling():
    # the canonical table keys on canonical labels; lookup must go via variants
    assert trait_class("age*sex") == "covariate"


def test_unknown_trait_is_not_guessed():
    assert trait_class("some trait nobody catalogued") is None


def test_finding_carries_its_trait_class():
    f = EwasCatalogProvider()._finding("cg1", _row("Tissue"))
    assert f.detail["trait_class"] == "covariate"


def test_finding_without_a_known_class_omits_the_key():
    f = EwasCatalogProvider()._finding("cg1", _row("some trait nobody catalogued"))
    assert "trait_class" not in f.detail
