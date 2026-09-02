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


_SHORT = {"age": "Epigenetic age", "bmi": "Body mass index", "crp": "C-reactive protein (inflammation marker)",
          "smoking": "Smoking", "sex": "Sex", "copd": "COPD", "eosinophilia": "Eosinophil count",
          "atopy": "Allergic sensitisation", "hdl_cholesterol": "HDL cholesterol",
          "total_cholesterol": "Total cholesterol", "body_fat": "Body fat", "waist_hip_ratio": "Waist-hip ratio",
          "gestational_age": "Gestational age at birth", "alcohol_consumption": "Alcohol consumption",
          "type_2_diabetes": "Type 2 diabetes", "chronic_pain": "Chronic pain", "hiv_infection": "HIV infection",
          "rheumatoid_arthritis": "Rheumatoid arthritis", "schizophrenia": "Schizophrenia"}


def short_label(d: dict, copy: dict) -> str:
    """A label a person can scan: the protein's name for a protein-level row,
    a short form for a known trait, else the trait with any parenthetical cut."""
    key = d.get("copy_key") or ""
    if key == "_protein_level" or d.get("protein"):
        subject = str(d.get("subject") or "")
        name = subject.replace("blood level of protein", "").strip() or str(d.get("protein") or "a protein")
        return f"{name} (a blood protein)"
    if key in _SHORT:
        return _SHORT[key]
    label = str(copy.get("label") or d.get("trait") or "this trait")
    return label.split(" (")[0].strip() or label


def _num(v):
    try:
        return None if v in (None, "") else float(v)
    except (TypeError, ValueError):
        return None


def _ewas(f: Finding) -> Interpretation:
    d = f.detail or {}
    copy = _copy_table().get(d.get("copy_key") or "", {}) or {}
    label = short_label(d, copy)
    d["short_label"] = label
    # The one clause worth a sentence: which way methylation moves with the
    # trait, when the studies say. Counts, tissues and p-values become chips.
    direction = d.get("direction")
    beta = _num(d.get("beta"))
    if direction == "mixed":
        clause = f"methylation here differs with {label.lower() if label[:1].isupper() and ' (' not in label else label}, but the studies disagree on which way"
    elif beta is not None and beta > 0:
        clause = "methylation here rises with it"
    elif beta is not None and beta < 0:
        clause = "methylation here falls with it"
    else:
        clause = "methylation here differs with it"
    found = f"{label} — {clause}."
    can = ("Group patterns at this site. Not a measurement of you and not a prediction.")
    cites = [ChainLink(kind="paper", label=f"PMID {pm}", id=f"PMID:{pm}",
                       url=f"https://pubmed.ncbi.nlm.nih.gov/{pm}/") for pm in (f.pmids or [])]
    return Interpretation(found=found, can_mean=can, how_sure="", next_step="",
                          condition=None, condition_ids=[], zygosity=None, citations=cites,
                          copy_version="trait_copy" if copy else "inline", reviewed_by=[])


def _gdc(f: Finding) -> Interpretation:
    d = f.detail or {}
    proj = str(d.get("project") or "a TCGA project")
    delta = _num(d.get("delta_beta"))
    d["short_label"] = proj.replace("TCGA-", "")
    found = f"Differs in {proj} tumour tissue" + (f" ({delta:+.2f})." if delta is not None else ".")
    can = "Tumour tissue from other people. Not a cancer test."
    url = f.link or f"https://portal.gdc.cancer.gov/projects/{proj}"
    return Interpretation(found=found, can_mean=can, how_sure="", next_step="",
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
