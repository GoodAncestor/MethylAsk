# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Compose what a methylation finding means, from its fields and the trait copy.

EWAS rows are population associations; TCGA rows are tumour-versus-normal
summaries. Both are group biology. The composer says so in each part, and
never turns the reader's own reading into a verdict.
"""
from __future__ import annotations
from biocore.providers.base import Finding, Interpretation, ChainLink
from .traits import _copy_table


def _ewas(f: Finding) -> Interpretation:
    d = f.detail or {}
    copy = _copy_table().get(d.get("copy_key") or "", {}) or {}
    trait = copy.get("label") or d.get("trait") or "this trait"
    reading = d.get("your reading")
    found = f"Studies link methylation at {f.marker} to {trait}."
    if reading is not None:
        found += f" Your reading at this site is {float(reading):.2f} on a 0 to 1 scale."
    can = copy.get("what_it_is_not") or (
        "This is a pattern seen across groups of people in research. "
        "It is not a measurement of your trait and not a prediction.")
    how = []
    p = d.get("p")
    if p not in (None, ""):
        try:
            how.append(f"The association reached p = {float(p):.0e}.")
        except (TypeError, ValueError):
            pass
    n = d.get("n")
    if n not in (None, ""):
        try:
            how.append(f"The study read {int(float(n)):,} people.")
        except (TypeError, ValueError):
            pass
    if d.get("tissue"):
        how.append(f"Tissue: {d['tissue']}.")
    if copy.get("typical_evidence"):
        how.append(str(copy["typical_evidence"]))
    cites = [ChainLink(kind="paper", label=f"PMID {pm}", id=f"PMID:{pm}",
                       url=f"https://pubmed.ncbi.nlm.nih.gov/{pm}/") for pm in (f.pmids or [])]
    return Interpretation(found=found, can_mean=str(can), how_sure=" ".join(how), next_step="",
                          condition=None, condition_ids=[], zygosity=None, citations=cites,
                          copy_version="trait_copy" if copy else "inline", reviewed_by=[])


def _gdc(f: Finding) -> Interpretation:
    d = f.detail or {}
    proj = d.get("project") or "a TCGA project"
    delta = d.get("delta_beta")
    found = f"In {proj}, methylation at {f.marker} differed between tumour tissue and normal tissue"
    found += f" by {float(delta):+.2f}." if delta is not None else "."
    can = ("This describes tumour tissue from other people. It is not a cancer test, "
           "and your reading was not compared with any diagnostic threshold.")
    nt, nn = int(d.get("n_tumor") or 0), int(d.get("n_normal") or 0)
    how = f"Tumour samples: {nt:,}. Normal samples: {nn:,}. Summary values, not individual cases."
    if d.get("sampling"):
        how += f" {d['sampling']}"
    url = f.link or f"https://portal.gdc.cancer.gov/projects/{proj}"
    return Interpretation(found=found, can_mean=can, how_sure=how, next_step="",
                          condition=None, condition_ids=[], zygosity=None,
                          citations=[ChainLink(kind="assertion", label=f"GDC {proj}", url=url)],
                          copy_version="inline", reviewed_by=[])


def interpret(findings: list[Finding]) -> int:
    """Fill interpretation and evidence_chain in place. Returns how many were filled."""
    n = 0
    for f in findings:
        if f.source == "ewas_catalog":
            f.interpretation = _ewas(f)
        elif f.source == "gdc":
            f.interpretation = _gdc(f)
        else:
            continue
        chain = [ChainLink(kind="variant", label=f.marker, url=f.link)]
        if (f.detail or {}).get("trait"):
            chain.append(ChainLink(kind="trait", label=str(f.detail["trait"]),
                                   id=f.detail.get("efo"), url=f.detail.get("efo")))
        f.evidence_chain = chain + list(f.interpretation.citations)
        n += 1
    return n
