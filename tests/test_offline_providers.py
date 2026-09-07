"""Native whole-genome analysis must not fall back to per-probe web requests."""
import sqlite3

from biocore.providers.base import Health
from methylask.providers import ewas_mirror
from methylask.providers.ewas_catalog import EwasCatalogProvider
from methylask.providers.clinvar import ClinVarProvider
from methylask.providers.gdc import GdcProvider


def test_offline_provider_status_and_lookups_never_network(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("offline provider reached network")
    monkeypatch.setattr("urllib.request.urlopen", forbidden)
    monkeypatch.setattr(ewas_mirror, "MIRROR_DB", tmp_path / "absent.db")
    providers = [EwasCatalogProvider(offline_only=True), ClinVarProvider(offline_only=True),
                 GdcProvider(summary_path=str(tmp_path / "absent-gdc.db"), offline_only=True)]
    for provider in providers:
        assert provider.status().health == Health.UNAVAILABLE
        assert provider.get_many(["cg00000001", "rs328"]) == {"cg00000001": [], "rs328": []}


def test_batched_ewas_uses_real_local_database(tmp_path, monkeypatch):
    path = tmp_path / "ewas.db"
    with sqlite3.connect(path) as con:
        con.execute("CREATE TABLE findings (cpg TEXT, trait TEXT, gene TEXT, beta TEXT, se TEXT, "
                    "p TEXT, n TEXT, tissue TEXT, methylation_array TEXT, chrpos TEXT, pmid TEXT, efo TEXT)")
        con.execute("INSERT INTO findings VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    ("cg00000001", "age", "TEST", "0.1", "0.01", "1e-10", "2000", "blood", "450K", "1:10", "123", ""))
    monkeypatch.setattr(ewas_mirror, "MIRROR_DB", path)
    provider = EwasCatalogProvider(offline_only=True)
    assert provider.status().health == Health.OK
    markers = [f"cg{i:08d}" for i in range(1100)]
    results = provider.get_many(markers)
    assert len(results) == 1100
    assert len(results["cg00000001"]) == 1
    assert results["cg00000002"] == []
    assert results["cg00000001"][0].detail["tissue"] == "blood"
