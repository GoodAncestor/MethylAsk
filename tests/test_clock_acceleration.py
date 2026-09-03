import csv
import math

import pytest

import methylask.clocks as clocks


def _clock_file(tmp_path, name, intercept, weights):
    path = tmp_path / f"{name}.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Probe", "Coefficient"])
        writer.writerow(["Intercept", intercept])
        writer.writerows(weights.items())
    return path


def test_linear_clock_acceleration_and_signed_contributions(tmp_path, monkeypatch):
    _clock_file(
        tmp_path,
        "Hannum2013_Blood",
        10,
        {"cg1": 2, "cg2": -1, "cg3": 0.5},
    )
    monkeypatch.setattr(clocks, "_CLOCK_DIR", tmp_path)

    result = clocks.Clock("Hannum2013_Blood").predict(
        {"cg1": 0.5, "cg2": 0.2, "cg3": 0.4},
        tissue="blood",
        age=8,
    )

    assert result.age == pytest.approx(11.0)
    assert result.acceleration == pytest.approx(3.0)
    assert result.contributions == [
        ("cg1", 2.0, 0.5, 1.0, 1.0),
        ("cg2", -1.0, 0.2, -0.2, -0.2),
        ("cg3", 0.5, 0.4, 0.2, 0.2),
    ]
    assert clocks.top_contributions(result, n=2) == result.contributions[:2]


@pytest.mark.parametrize(
    ("intercept", "beta", "expected_slope"),
    [
        (-1.0, 0.5, 21 * math.exp(-0.5)),
        (0.1, 0.5, 21.0),
    ],
)
def test_horvath_contribution_uses_the_local_inverse_slope(
    tmp_path, monkeypatch, intercept, beta, expected_slope
):
    _clock_file(
        tmp_path,
        "Horvath2013_PanTissue",
        intercept,
        {"cg1": 1, "cg2": 0, "cg3": 0},
    )
    monkeypatch.setattr(clocks, "_CLOCK_DIR", tmp_path)

    result = clocks.Clock("Horvath2013_PanTissue").predict(
        {"cg1": beta, "cg2": 0.0, "cg3": 0.0}, age=30
    )

    assert result.contributions[0][4] == pytest.approx(beta * expected_slope)
    assert result.acceleration == pytest.approx(result.age - 30)


def test_acceleration_requires_an_age_and_a_valid_clock(tmp_path, monkeypatch):
    _clock_file(tmp_path, "Hannum2013_Blood", 10, {"cg1": 2, "cg2": 1, "cg3": 1})
    monkeypatch.setattr(clocks, "_CLOCK_DIR", tmp_path)
    clock = clocks.Clock("Hannum2013_Blood")

    assert clock.predict({"cg1": 0.5}, tissue="blood").acceleration is None
    invalid = clock.predict({"cg1": 0.5}, tissue="buccal", age=8)
    assert invalid.valid is False
    assert invalid.acceleration is None


def test_run_all_passes_age_to_each_clock(monkeypatch):
    calls = []

    class FakeClock:
        def __init__(self, name, min_coverage):
            self.name = name

        def predict(self, betas, tissue=None, age=None):
            calls.append((self.name, age))
            return self.name

    monkeypatch.setattr(clocks, "Clock", FakeClock)

    clocks.run_all({}, age=52)

    assert calls and all(age == 52 for _name, age in calls)
