"""methylask CLI: status | refresh | report (docs/DESIGN.md §6)."""
from __future__ import annotations
import argparse, sys
from .providers.registry import Registry
from .providers.ewas_catalog import EwasCatalogProvider


def _default_registry() -> Registry:
    return Registry([EwasCatalogProvider()])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="methylask")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status", help="provider health + cache ages")
    pr = sub.add_parser("refresh", help="build/refresh local caches")
    pr.add_argument("--provider", default="all")
    rp = sub.add_parser("report", help="produce a report from a sample file")
    rp.add_argument("sample")
    rp.add_argument("--pdf", action="store_true")
    args = ap.parse_args(argv)

    reg = _default_registry()
    if args.cmd == "status":
        for s in reg.status():
            print(f"{s.name:20s} {s.health.value:12s} {s.note or ''}")
        return 0
    if args.cmd == "refresh":
        print("refresh: not yet implemented in scaffold")
        return 0
    if args.cmd == "report":
        from .ingest.beta_matrix import read_beta_matrix
        sample = read_beta_matrix(args.sample)
        rep = reg.annotate(sample.markers)
        print(f"{len(sample.markers)} markers -> {len(rep.all_findings())} findings")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
