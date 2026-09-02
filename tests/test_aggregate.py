from biocore.providers.base import Category, Finding, Tier

from methylask.aggregate import aggregate_by_trait


def _row(
    pmid,
    beta,
    p,
    n,
    tissue="whole blood",
    trait="Body mass index",
    marker="cg1",
    tier=Tier.MODERATE,
):
    return Finding(
        marker=marker,
        source="ewas_catalog",
        description="x",
        tier=tier,
        categories=[Category.CLINICAL],
        detail={
            "trait": trait,
            "copy_key": "bmi",
            "beta": beta,
            "p": p,
            "n": n,
            "tissue": tissue,
            "your reading": 0.41,
            "topic": "metabolic",
        },
        pmids=[pmid],
    )


def test_rows_for_one_marker_and_trait_become_one_finding():
    findings = [
        _row("1", 0.02, 1e-12, 5387, tier=Tier.ROBUST),
        _row("2", 0.015, 1e-6, 1200),
        _row("3", 0.03, 1e-9, 3000, tissue="cord blood"),
    ]
    out = aggregate_by_trait(findings, "blood")
    assert len(out) == 1
    finding = out[0]
    assert finding.tier == Tier.ROBUST
    assert sorted(finding.pmids) == ["1", "2", "3"]
    detail = finding.detail
    assert detail["n_studies"] == 3
    assert detail["n_participants"] == 9587
    assert detail["direction"] == "consistent"
    assert detail["p"] == 1e-12
    assert detail["n"] == 9587
    assert detail["your reading"] == 0.41
    assert sorted(detail["tissues"]) == ["cord blood", "whole blood"]
    assert detail["tissue_supported"] is True
    assert len(detail["rows"]) == 3


def test_different_traits_and_markers_stay_apart_and_mixed_direction_is_named():
    findings = [
        _row("1", 0.02, 1e-9, 100),
        _row("2", -0.02, 1e-9, 100),
        _row("3", 0.02, 1e-9, 100, trait="Age"),
        _row("4", 0.02, 1e-9, 100, marker="cg2"),
    ]
    out = aggregate_by_trait(findings, None)
    assert len(out) == 3
    bmi = next(
        finding
        for finding in out
        if finding.detail["trait"] == "Body mass index" and finding.marker == "cg1"
    )
    assert bmi.detail["direction"] == "mixed"
    assert bmi.detail["n_studies"] == 2


def test_non_ewas_findings_pass_through_in_place():
    other = Finding(
        marker="cg1",
        source="gdc",
        description="x",
        tier=Tier.ROBUST,
        categories=[Category.REFERENCE],
    )
    out = aggregate_by_trait([_row("1", 0.02, 1e-9, 100), other], None)
    assert other in out
    assert len(out) == 2
