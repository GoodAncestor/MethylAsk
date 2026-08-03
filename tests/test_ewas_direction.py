"""EWAS finding descriptions state the direction of the association.

"linked to age" tells a reader nothing about how it is linked. The sign of the
published effect size is the direction, and it was already being carried in
detail as a bare signed number behind a disclosure.
"""
from methylask.providers.ewas_catalog import EwasCatalogProvider


def _row(beta, trait="age", **kw):
    r = {"trait": trait, "gene": "ERGIC3", "beta": beta, "se": 0.0002,
         "p": 1e-30, "n": 2338, "tissue": "Whole blood",
         "methylation_array": "450k", "chrpos": "20:35", "pmid": "33450751"}
    r.update(kw)
    return r


def test_negative_effect_reads_as_lower_methylation():
    f = EwasCatalogProvider()._finding("cg00017842", _row(-0.003))
    assert "lower methylation" in f.description
    assert "age" in f.description


def test_positive_effect_reads_as_higher_methylation():
    f = EwasCatalogProvider()._finding("cg00017842", _row(0.004))
    assert "higher methylation" in f.description


def test_missing_effect_size_falls_back_to_the_plain_wording():
    # no published effect -> no direction claimed
    f = EwasCatalogProvider()._finding("cg00017842", _row(None))
    assert f.description == "linked to age"


def test_effect_size_and_sample_size_are_kept_for_display():
    f = EwasCatalogProvider()._finding("cg00017842", _row(-0.003))
    assert f.detail["beta"] == -0.003
    assert f.detail["n"] == 2338
