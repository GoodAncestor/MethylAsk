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


# --- the copy key must survive humanization --------------------------------

def test_protein_finding_keeps_a_usable_copy_key():
    """humanize_trait rewrites 'P01023' to 'Alpha-2-macroglobulin' before it is
    stored, so a renderer looking up copy from detail['trait'] finds nothing for
    any protein — 68.7% of the catalog. The key is resolved from the RAW trait
    at provider time and carried on the finding."""
    from methylask.providers.ewas_catalog import EwasCatalogProvider
    from methylask.traits import trait_copy_key, _copy_table
    f = EwasCatalogProvider()._finding("cg1", {
        "trait": "P01023", "beta": -0.01, "n": 2000, "p": 1e-20, "se": 0.001,
        "tissue": "Whole blood", "methylation_array": "450k", "chrpos": "1:1",
        "pmid": "1", "gene": "A2M"})
    assert f.detail["trait"] == "Alpha-2-macroglobulin"      # display name kept
    assert f.detail["copy_key"] == "_protein_level"          # copy still reachable
    assert _copy_table()[f.detail["copy_key"]]["label"] == "Blood level of a protein"


def test_named_trait_finding_carries_its_copy_key():
    from methylask.providers.ewas_catalog import EwasCatalogProvider
    f = EwasCatalogProvider()._finding("cg1", {
        "trait": "BMI", "beta": 0.01, "n": 2000, "p": 1e-20, "se": 0.001,
        "tissue": "Whole blood", "methylation_array": "450k", "chrpos": "1:1",
        "pmid": "1", "gene": "X"})
    assert f.detail["copy_key"] == "bmi"


def test_uncurated_trait_carries_no_copy_key():
    from methylask.providers.ewas_catalog import EwasCatalogProvider
    f = EwasCatalogProvider()._finding("cg1", {
        "trait": "Fractional exhaled nitric oxide (FeNO)", "beta": 0.01, "n": 20,
        "p": 0.01, "se": 0.001, "tissue": "Whole blood", "methylation_array": "450k",
        "chrpos": "1:1", "pmid": "1", "gene": "X"})
    assert "copy_key" not in f.detail


def test_trait_copy_key_is_none_for_unknown():
    from methylask.traits import trait_copy_key
    assert trait_copy_key("nothing at all") is None
