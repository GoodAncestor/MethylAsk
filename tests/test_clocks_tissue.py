"""Pin the clock tissue-mismatch + plausibility flagging.

A blood-trained clock applied to a non-blood tissue must be flagged (valid=False)
rather than presenting its number; a pan-tissue clock is never tissue-flagged;
an implausible age is flagged. Uses each clock's own CpG set so coverage is full
and the only variable under test is the tissue/plausibility logic.
"""
import csv
from pathlib import Path
import methylask.clocks as ck

CLOCK_DIR = Path(ck.__file__).parent / "data" / "reference" / "clocks"


def _betas_for(clock_name, value=0.5):
    """Full-coverage betas for one clock's CpGs (excludes Intercept)."""
    betas = {}
    with open(CLOCK_DIR / f"{clock_name}.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            if row["Probe"] != "Intercept":
                betas[row["Probe"]] = value
    return betas


def test_pan_tissue_never_tissue_flagged():
    betas = _betas_for("Horvath2013_PanTissue")
    r = ck.Clock("Horvath2013_PanTissue").predict(betas, tissue="buccal")
    assert r.trained_tissue == "pan"
    assert r.tissue_mismatch is False   # pan clocks apply to any tissue


def test_blood_clock_flagged_on_buccal():
    betas = _betas_for("Hannum2013_Blood")
    r = ck.Clock("Hannum2013_Blood").predict(betas, tissue="buccal")
    assert r.tissue_mismatch is True
    assert r.valid is False
    assert "not valid for this sample type" in r.note


def test_blood_clock_ok_on_blood():
    betas = _betas_for("Hannum2013_Blood")
    r = ck.Clock("Hannum2013_Blood").predict(betas, tissue="blood")
    assert r.tissue_mismatch is False


def test_no_tissue_arg_never_mismatches():
    # backward compatible: omitting tissue must not flag anything as mismatch
    betas = _betas_for("Levine2018_PhenoAge")
    r = ck.Clock("Levine2018_PhenoAge").predict(betas)
    assert r.tissue_mismatch is False


def test_implausible_age_flagged():
    # force an implausible predictor by giving PhenoAge extreme betas
    betas = _betas_for("Levine2018_PhenoAge", value=100.0)
    r = ck.Clock("Levine2018_PhenoAge").predict(betas, tissue="blood")
    # extreme input -> age outside human-plausible range -> flagged, valid False
    assert r.implausible is True
    assert r.valid is False


def test_run_all_threads_tissue():
    betas = {}
    for name in ck.available():
        betas.update(_betas_for(name))
    results = ck.run_all(betas, tissue="saliva")
    by = {r.clock: r for r in results}
    assert by["Horvath2013_PanTissue"].tissue_mismatch is False
    assert by["Hannum2013_Blood"].tissue_mismatch is True
