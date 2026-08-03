"""Plain-language copy lookup for a trait, keyed from the raw catalog string."""
from methylask.traits import trait_copy

FIELDS = ("what_it_is", "what_an_association_means", "what_it_is_not", "typical_evidence")


def test_raw_variant_resolves_to_its_copy():
    c = trait_copy("BMI")
    assert c["label"]
    assert all(c.get(f) for f in FIELDS)


def test_alias_forms_of_one_concept_share_copy():
    # "Incident COPD" and "Prevalent COPD (Self-report)" are separate concepts for
    # counting, but the explanation of what COPD IS is common to both
    assert trait_copy("Incident COPD")["label"] == trait_copy("Prevalent COPD (Self-report)")["label"]


def test_protein_traits_get_the_class_wide_entry():
    # 68.7% of associations are protein levels; they share one explanation
    c = trait_copy("P02748")
    assert c["label"] == "Blood level of a protein"
    assert "pQTM" in c["what_an_association_means"] or "protein" in c["what_an_association_means"]


def test_written_form_protein_also_gets_the_class_entry():
    assert trait_copy("PAPPA protein levels (SeqId = 4148-49)")["label"] == "Blood level of a protein"


def test_uncurated_trait_returns_nothing_rather_than_guessing():
    assert trait_copy("Fractional exhaled nitric oxide (FeNO)") == {}


def test_every_copy_entry_carries_pmids():
    from methylask.traits import _copy_table
    for k, v in _copy_table().items():
        assert v.get("pmids"), f"{k} has no pmids"
