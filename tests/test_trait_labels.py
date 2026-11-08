"""Outcome labels a reader sees never read as catalogue text.

The EWAS Catalog's trait strings are study spellings: "age", "age*sex",
"Age 4 vs age 0", "M32412_butyrylcarnitine". A report groups findings into
outcomes titled by the finding's trait label, so any of those spellings that
reaches the label is shown to a lay reader as the name of an outcome.
"""
import json
import re
from pathlib import Path

import pytest

from methylask.meaning import interpret
from methylask.providers.ewas_catalog import EwasCatalogProvider
from methylask.traits import display_label, trait_class, trait_copy_key

_TABLE = json.loads((Path(__file__).parent.parent / "methylask" / "data" / "reference"
                     / "trait_canonical.json").read_text())

# What catalogue text looks like. A rendered outcome label must match none.
CATALOGUE_PATTERNS = {
    # "age*sex"; a trailing "*" on a Metabolon compound name is not an operator
    "covariate operator": re.compile(r"\S\*\S"),
    "interaction term": re.compile(r"\binteraction\b", re.I),
    "age contrast": re.compile(r"\bage \d+\s+vs\.?\s+age \d+\b", re.I),
    "lower-case raw spelling": re.compile(r"^[a-z]{2,}(?=[\s,'-]|$)"),
    "Metabolon code": re.compile(r"^M\d{5}_"),
}

# Raw strings seen in the HG002 report's By outcome view (2026-09-06 render)
# or in the mirror, beyond the curated table's own variants.
_RAW_SEEN = ["age", "sex", "schizophrenia", "body mass index", "breast cancer",
             "fasting glucose", "fasting insulin", "psychosis", "monozygotic twinning",
             "alzheimers disease braak stage", "hepatic fat", "waist circumference-to-hip ratio",
             "homocysteine (MTHFR 677C>T)", "M32412_butyrylcarnitine",
             "M32654_3-dehydrocarnitine*", "Age 4 vs age 0", "Ageing"]


def _matches(label):
    return [name for name, pat in CATALOGUE_PATTERNS.items() if pat.search(label)]


def _row(trait):
    return {"trait": trait, "gene": "X", "beta": 0.01, "se": 0.001, "p": 1e-9, "n": 100,
            "tissue": "Whole blood", "methylation_array": "450k", "chrpos": "1:1", "pmid": "1"}


def _shown():
    """Every raw string a report would show (covariates are withheld)."""
    raws = [v for c in _TABLE if c["class"] != "covariate" for v in c["variants"]]
    return raws + _RAW_SEEN


@pytest.mark.parametrize("raw", _shown())
def test_no_rendered_label_reads_as_catalogue_text(raw):
    f = EwasCatalogProvider()._finding("cg1", _row(raw))
    interpret([f])
    rendered = [f.detail["trait"], f.detail.get("short_label", ""),
                f.interpretation.found.split(" — ")[0]]
    for label in rendered:
        assert not _matches(label), f"{raw!r} renders as {label!r}: {_matches(label)}"


def test_no_canonical_label_of_a_shown_concept_is_catalogue_text():
    bad = {c["canonical_label"]: _matches(c["canonical_label"])
           for c in _TABLE if c["class"] != "covariate" and _matches(c["canonical_label"])}
    assert not bad


@pytest.mark.parametrize("raw", ["age*sex", "Age-by-sex interaction",
                                 "Interaction between gestational age and sex",
                                 "sex-Braak stage interaction", "Tissue", "Fetal vs adult liver"])
def test_study_design_terms_are_covariates(raw):
    assert trait_class(raw) == "covariate"


def test_every_interaction_or_operator_term_in_the_table_is_a_covariate():
    for c in _TABLE:
        for v in c["variants"]:
            if CATALOGUE_PATTERNS["covariate operator"].search(v):
                assert c["class"] == "covariate", v
            if CATALOGUE_PATTERNS["interaction term"].search(v):
                assert c["class"] == "covariate", v


@pytest.mark.parametrize("raw", ["age", "Age"])
def test_age_spellings_share_one_age_outcome(raw):
    assert display_label(raw) == "Age"
    assert trait_copy_key(raw) == "age"


def test_age_contrast_reads_as_a_comparison_and_keeps_its_own_outcome():
    # cord blood at birth vs whole blood at age 4 (PMID 28056824): a specific
    # developmental comparison, not general ageing, so it neither merges into
    # Age nor borrows the age explanation
    assert display_label("Age 4 vs age 0") == "Age comparison: 4 years vs birth"
    assert trait_copy_key("Age 4 vs age 0") is None


def test_age_adjacent_concepts_do_not_fold_into_age():
    assert display_label("Ageing") == "Ageing"
    assert trait_copy_key("Ageing") is None
    assert display_label("Gestational age") == "Gestational age"
    assert trait_copy_key("Gestational age") == "gestational_age"
    assert trait_copy_key("Age at menarche") is None


def test_raw_lower_case_labels_are_capitalised():
    assert display_label("hepatic fat") == "Hepatic fat"
    assert display_label("schizophrenia") == "Schizophrenia"
    assert display_label("alzheimers disease braak stage") == "Alzheimer's disease Braak stage"


@pytest.mark.parametrize("name", ["miR-100-5p", "eGFR decline", "p53 status", "mtDNA copy number"])
def test_case_sensitive_names_are_left_alone(name):
    assert display_label(name) == name


def test_metabolon_codes_are_spelled_out():
    assert display_label("M32412_butyrylcarnitine") == \
        "Metabolite level: butyrylcarnitine (Metabolon ID M32412)"
    assert display_label("M32654_3-dehydrocarnitine*") == \
        "Metabolite level: 3-dehydrocarnitine* (Metabolon ID M32654)"
    # matches the curated table's own spelling of the same kind of readout
    assert display_label("M00054_tryptophan") == "Metabolite level: tryptophan (Metabolon ID M00054)"


def test_spelling_variants_aggregate_at_the_same_site():
    from methylask.aggregate import aggregate_by_trait
    p = EwasCatalogProvider()
    rows = [{**_row(r), "pmid": pmid} for r, pmid in (("age", "1"), ("Age", "2"), ("Age 4 vs age 0", "3"))]
    out = aggregate_by_trait([p._finding("cg1", r) for r in rows], None)
    by = {f.detail["trait"]: f for f in out}
    assert set(by) == {"Age", "Age comparison: 4 years vs birth"}
    assert by["Age"].detail["n_studies"] == 2
