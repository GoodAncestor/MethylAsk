"""Published reference values for a marker.

These are absolute population levels quoted from the literature (e.g. the mean
beta at cg05575921 in never-smokers), each carrying the PMID it came from. They
are NOT derivable from the EWAS mirror, which stores effect sizes — differences
between groups — not population levels. So this table is curated, and a marker
with no published reference must say so rather than be silently dropped.
"""
import json
import pytest
from methylask.reference import (load_reference_table, reference_for,
                                 describe_position, positions_for_sample)


@pytest.fixture
def table(tmp_path):
    p = tmp_path / "marker_reference.json"
    p.write_text(json.dumps({
        "cg05575921": {
            "gene": "AHRR",
            "references": [
                {"group": "never-smoker", "beta": 0.85, "sd": 0.03,
                 "tissue": "whole blood", "array": "450k", "n": 1793,
                 "pmid": "00000000"},
            ],
        },
    }))
    return load_reference_table(p)


def test_marker_with_a_published_reference_returns_it(table):
    refs = reference_for(table, "cg05575921")
    assert len(refs) == 1
    assert refs[0]["group"] == "never-smoker"
    assert refs[0]["pmid"] == "00000000"


def test_marker_without_a_published_reference_returns_empty(table):
    # not an error — absence of a reference value is itself reportable
    assert reference_for(table, "cg00000000") == []


def test_position_states_the_gap_and_cites_the_source(table):
    ref = reference_for(table, "cg05575921")[0]
    pos = describe_position(sample_beta=0.71, ref=ref)
    assert pos.reference_group == "never-smoker"
    assert pos.reference_beta == 0.85
    assert pos.pmid == "00000000"
    assert pos.delta == pytest.approx(-0.14)


def test_position_carries_the_statistic_type_so_copy_cannot_say_average(table):
    # most published values are medians or adjusted means; calling a median an
    # "average" in report copy misdescribes the source
    ref = dict(reference_for(table, "cg05575921")[0], stat="median")
    pos = describe_position(sample_beta=0.71, ref=ref)
    assert pos.stat == "median"


def test_position_carries_tissue_and_n_for_caveating(table):
    # a cord-blood or n=16 reference must be surfaceable at the point of display
    ref = reference_for(table, "cg05575921")[0]
    pos = describe_position(sample_beta=0.71, ref=ref)
    assert pos.tissue == "whole blood"
    assert pos.n == 1793


def test_beta_outside_zero_to_one_is_rejected_on_load(tmp_path):
    # a value like 85 is a percentage and -0.17 is an effect size; either pasted
    # into a beta field is a silent, large error
    p = tmp_path / "bad.json"
    p.write_text(json.dumps({
        "cg00000001": {"gene": "X", "references": [
            {"group": "cohort", "beta": 85.0, "tissue": "whole blood",
             "array": "450k", "n": 100, "pmid": "1"}]}}))
    with pytest.raises(ValueError, match="cg00000001"):
        load_reference_table(p)


def test_position_without_a_published_sd_has_no_sigma(table):
    # a sigma figure invented without a published SD would be a fabricated statistic
    ref = dict(reference_for(table, "cg05575921")[0])
    ref.pop("sd")
    pos = describe_position(sample_beta=0.71, ref=ref)
    assert pos.sigma is None


# --- whole-sample positioning, incl. the tissue gate -------------------------

def test_blood_sample_gets_positions_against_blood_references(table):
    out = positions_for_sample({"cg05575921": 0.93}, tissue="blood", table=table)
    assert len(out) == 1
    assert out[0].probe == "cg05575921"
    assert out[0].suppressed_reason is None
    assert len(out[0].positions) == 1


def test_buccal_sample_suppresses_blood_only_references(table):
    # a pediatric buccal 0.56 sits nearest the whole-blood CURRENT-SMOKER median.
    # Showing that comparison would read a child's cheek swab as a smoker, so a
    # tissue-mismatched reference is withheld, matching clocks.py's precedent.
    out = positions_for_sample({"cg05575921": 0.56}, tissue="buccal", table=table)
    assert len(out) == 1
    assert out[0].positions == []
    assert "buccal" in out[0].suppressed_reason
    assert "whole blood" in out[0].suppressed_reason


def test_probes_with_no_curated_reference_are_skipped(table):
    out = positions_for_sample({"cg99999999": 0.5}, tissue="blood", table=table)
    assert out == []


def test_unknown_sample_tissue_does_not_suppress(table):
    # omitting tissue must not silently withhold everything
    out = positions_for_sample({"cg05575921": 0.93}, tissue=None, table=table)
    assert out[0].positions and out[0].suppressed_reason is None


# --- human-readable meaning --------------------------------------------------

def test_marker_meaning_exposes_label_and_the_two_caveat_fields():
    from methylask.reference import marker_meaning
    t = load_reference_table()            # the shipped, curated table
    m = marker_meaning(t, "cg05575921")
    assert m["label"] == "Tobacco smoke exposure"
    assert "AHRR" in m["what_was_read"]
    assert m["what_it_is_not"]


def test_marker_meaning_is_empty_for_an_uncurated_marker():
    from methylask.reference import marker_meaning
    assert marker_meaning({}, "cg99999999") == {}


def test_every_curated_marker_has_human_copy():
    # a marker with a reference value but no plain-language label would render
    # as a bare probe id and a bare number, which is what we are fixing
    from methylask.reference import marker_meaning
    t = load_reference_table()
    for probe in t:
        m = marker_meaning(t, probe)
        assert m.get("label"), f"{probe} has no label"
        assert m.get("what_was_read"), f"{probe} has no what_was_read"
        assert m.get("what_it_is_not"), f"{probe} has no what_it_is_not"
