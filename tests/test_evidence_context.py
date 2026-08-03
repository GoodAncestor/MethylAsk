"""Replication context derived from EWAS mirror rows.

Morgan's prototype hand-writes a replication line per card ("40+ cohorts,
meta-analysed" vs "3 cohorts, mixed direction"). Those are derivable from the
mirror's per-association rows, so the robustness statement is computed from the
evidence rather than assigned by a human.
"""
from methylask.evidence import summarize_replication


def _row(pmid, beta, n=1000, tissue="whole blood", trait="smoking"):
    return {"trait": trait, "gene": "AHRR", "beta": beta, "se": 0.01,
            "p": 1e-20, "n": n, "tissue": tissue,
            "methylation_array": "450k", "chrpos": "5:373378",
            "pmid": pmid, "efo": "EFO_0004318"}


def test_studies_agreeing_in_sign_are_consistent():
    rows = [_row("1", -0.08), _row("2", -0.06), _row("3", -0.09)]
    ctx = summarize_replication(rows)
    assert ctx.direction == "consistent"


def test_studies_disagreeing_in_sign_are_mixed():
    rows = [_row("1", -0.08), _row("2", 0.06), _row("3", -0.09)]
    ctx = summarize_replication(rows)
    assert ctx.direction == "mixed"


def test_one_study_reporting_several_traits_counts_once():
    # the same paper contributes a row per trait; that is one cohort, not three
    rows = [_row("1", -0.08, trait="smoking"),
            _row("1", -0.07, trait="cigarettes per day"),
            _row("2", -0.06, trait="smoking")]
    ctx = summarize_replication(rows)
    assert ctx.n_studies == 2


def test_participants_counted_once_per_study_not_per_row():
    # summing the n column blindly would double-count this study's 1000 people
    rows = [_row("1", -0.08, n=1000, trait="smoking"),
            _row("1", -0.07, n=1000, trait="cigarettes per day"),
            _row("2", -0.06, n=500, trait="smoking")]
    ctx = summarize_replication(rows)
    assert ctx.n_participants == 1500


def test_tissue_is_reported_as_the_set_of_tissues_studied():
    rows = [_row("1", -0.08, tissue="whole blood"),
            _row("2", -0.06, tissue="whole blood"),
            _row("3", -0.09, tissue="buccal")]
    ctx = summarize_replication(rows)
    assert ctx.tissues == ["buccal", "whole blood"]   # sorted, deduplicated


def test_sample_tissue_absent_from_the_literature_is_unsupported():
    rows = [_row("1", -0.08, tissue="whole blood")]
    ctx = summarize_replication(rows, sample_tissue="saliva")
    assert ctx.tissue_supported is False


def test_sample_tissue_present_in_the_literature_is_supported():
    rows = [_row("1", -0.08, tissue="whole blood")]
    ctx = summarize_replication(rows, sample_tissue="blood")
    assert ctx.tissue_supported is True
