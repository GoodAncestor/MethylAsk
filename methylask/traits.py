"""Trait classification + humanization for EWAS findings.

Two jobs, both turning raw EWAS trait strings into something a person can read
and filter:

1. classify_topic(trait) -> a subject tag (aging, metabolic, cardiovascular,
   immune, cancer, neuro, reproductive, lifestyle, proteomic, other) so the
   report can offer a "show me things about age" filter.

2. humanize_trait(trait) -> (label, kind, accession) — a bare UniProt accession
   like 'P02750' or 'A0A075B7B8' is recognised as a protein-level measurement and
   flagged for name resolution + linkout, rather than shown as a meaningless code.

Both are heuristic (keyword / pattern based): good, not perfect. The topic is a
navigation aid, never a clinical claim.
"""
from __future__ import annotations
import re, json, functools
from pathlib import Path

_UNIPROT_CACHE = Path(__file__).parent / "data" / "reference" / "uniprot_names.json"


@functools.lru_cache(maxsize=1)
def _uniprot_names() -> dict:
    """Bundled accession -> {name, gene} cache (resolved offline from UniProt).
    A production build refreshes this; unknown accessions fall back to the code."""
    try:
        return json.loads(_UNIPROT_CACHE.read_text())
    except Exception:
        return {}


def protein_name(accession: str) -> str | None:
    """Human protein name for a UniProt accession, or None if not in the cache."""
    rec = _uniprot_names().get((accession or "").strip())
    return rec.get("name") if rec else None

# UniProt accession: canonical (P12345 / Q9Y6K9) and the newer long form (A0A075B7B8)
_UNIPROT = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")

# topic keyword rules, checked in order (first match wins). Lowercased substring match.
_TOPICS = [
    ("aging",          ["age", "aging", "ageing", "gestational age", "longevity", "epigenetic age"]),
    ("cancer",         ["cancer", "tumor", "tumour", "carcinoma", "leukemia", "leukaemia",
                        "lymphoma", "melanoma", "neoplasm", "malignan"]),
    ("metabolic",      ["bmi", "body fat", "body mass", "obesit", "waist", "cholesterol", "hdl",
                        "ldl", "triglyceride", "glucose", "insulin", "diabet", "metaboli",
                        "lipid", "adiposity"]),
    ("cardiovascular", ["blood pressure", "hypertension", "cardiovascular", "coronary", "heart",
                        "cardiac", "stroke", "atheroscler"]),
    ("immune",         ["c-reactive protein", "crp", "inflammat", "immune", "cytokine",
                        "interleukin", "leukocyte", "autoimmun", "allerg", "asthma",
                        "rheumatoid", "colitis", "inflammatory bowel", "lupus", "hiv",
                        "eosinophil", "psoriasis", "arthritis"]),
    ("respiratory",    ["copd", "fev1", "fev-1", "lung function", "pulmonary", "respirat"]),
    ("neuro",          ["cognit", "alzheimer", "dementia", "parkinson", "depress", "schizophren",
                        "neuro", "brain", "psychiat"]),
    ("reproductive",   ["pregnan", "birth weight", "maternal", "fetal", "gestation", "menopaus",
                        "puberty", "fertil"]),
    ("lifestyle",      ["smoking", "smoke", "tobacco", "alcohol", "drinking", "diet", "exercise",
                        "physical activity", "sleep"]),
]


def classify_topic(trait: str) -> str:
    """Best-effort subject tag for a trait string."""
    t = (trait or "").lower().strip()
    if not t:
        return "other"
    # a bare accession, or a "<X> protein levels (SeqId=...)" string, is a
    # protein-level measurement (plasma proteomics EWAS)
    if _UNIPROT.match(trait.strip()) or "protein level" in t or "seqid" in t:
        return "proteomic"
    for topic, needles in _TOPICS:
        for n in needles:
            if n in t:
                return topic
    return "other"


def is_uniprot(trait: str) -> bool:
    return bool(_UNIPROT.match((trait or "").strip()))


def humanize_trait(trait: str) -> tuple[str, str, str | None]:
    """Return (display_label, kind, accession).

    kind is 'protein' for a UniProt accession (label stays the accession until a
    name is resolved by the caller), else 'trait' with the label unchanged.
    """
    t = (trait or "").strip()
    if is_uniprot(t):
        # resolve the accession to a readable protein name where we can; the
        # accession is still returned so the renderer can link out to UniProt
        name = protein_name(t)
        return (name or t, "protein", t)
    # light tidy: collapse whitespace, keep the study's wording otherwise
    return (re.sub(r"\s+", " ", t), "trait", None)


_CANONICAL = Path(__file__).parent / "data" / "reference" / "trait_canonical.json"


@functools.lru_cache(maxsize=1)
def _class_by_variant() -> dict:
    """raw trait string (lowercased) -> class, from the curated vocabulary.

    Keyed by VARIANT, not canonical label: the catalog's rows carry the raw
    spelling ("age*sex"), while the table names concepts ("Age x Sex interaction
    term"). Looking up by label would silently match nothing.
    """
    try:
        with open(_CANONICAL) as fh:
            table = json.load(fh)
    except (OSError, ValueError):
        return {}
    out = {}
    for c in table:
        for v in c.get("variants", []):
            out[str(v).strip().lower()] = c.get("class")
    return out


def trait_class(trait: str) -> str | None:
    """'health_trait' | 'covariate' | 'protein_level' | 'other', or None.

    None means the trait is not in the curated top-400 — the great majority of
    the 6,515 distinct strings. Absence is not evidence that a trait is safe to
    show; it only means nobody has classified it yet.
    """
    return _class_by_variant().get((trait or "").strip().lower())


_COPY = Path(__file__).parent / "data" / "reference" / "trait_copy.json"
_COPY_ALIAS = Path(__file__).parent / "data" / "reference" / "trait_copy_alias.json"


@functools.lru_cache(maxsize=1)
def _copy_table() -> dict:
    try:
        with open(_COPY) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


@functools.lru_cache(maxsize=1)
def _canonical_label_by_variant() -> dict:
    """raw trait string (lowercased) -> canonical_label."""
    try:
        with open(_CANONICAL) as fh:
            table = json.load(fh)
    except (OSError, ValueError):
        return {}
    return {str(v).strip().lower(): c["canonical_label"]
            for c in table for v in c.get("variants", [])}


@functools.lru_cache(maxsize=1)
def _copy_key_by_label() -> dict:
    """canonical_label -> copy key, via the curated alias map.

    Explicit rather than fuzzy: normalising labels to keys silently failed on 8 of
    22 concepts — including BMI, CRP, smoking, COPD and type 2 diabetes, i.e. most
    of the volume. A near-miss here attaches the wrong explanation to a trait, so
    the mapping is curated and asserted, never guessed.
    """
    try:
        with open(_COPY_ALIAS) as fh:
            alias = json.load(fh)
    except (OSError, ValueError):
        return {}
    return {label: key for key, labels in alias.items() for label in labels}


def trait_copy_key(trait: str) -> str | None:
    """Copy-table key for a RAW catalog trait string, or None.

    Must be resolved from the raw string: humanize_trait rewrites "P01023" to
    "Alpha-2-macroglobulin" before the trait is stored on a finding, and the
    rewritten name resolves to nothing — which would silently strip copy from
    every protein trait, i.e. 68.7% of the catalog.
    """
    t = (trait or "").strip()
    if not t:
        return None
    if is_uniprot(t) or re.search(r"protein levels?\s*\(SeqId", t, re.I):
        return "_protein_level" if "_protein_level" in _copy_table() else None
    label = _canonical_label_by_variant().get(t.lower())
    if not label:
        return None
    key = _copy_key_by_label().get(label)
    return key if key in _copy_table() else None


def trait_copy(trait: str) -> dict:
    """Plain-language copy for a raw trait string, or {} when none is curated."""
    return _copy_table().get(trait_copy_key(trait) or "", {})
