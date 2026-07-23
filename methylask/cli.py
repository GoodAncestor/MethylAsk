"""methylask CLI: status | refresh | report (docs/DESIGN.md §6)."""
from __future__ import annotations
import argparse, sys
from biocore.providers.registry import Registry
from .providers.ewas_catalog import EwasCatalogProvider
from .providers.clinvar import ClinVarProvider
from .providers.gdc import GdcProvider


def _default_registry() -> Registry:
    return Registry([
        EwasCatalogProvider(),
        ClinVarProvider(),
        GdcProvider(),
    ])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="methylask")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status", help="provider health + cache ages")
    pr = sub.add_parser("refresh", help="build/refresh local caches")
    pr.add_argument("--provider", default="all")
    rp = sub.add_parser("report", help="produce a report from a sample file")
    rp.add_argument("sample")
    rp.add_argument("--out", default="report.html")
    rp.add_argument("--pdf", action="store_true")
    rp.add_argument("--max-markers", type=int, default=200,
                    help="cap on markers annotated via live API (guards against "
                         "full-array hangs until the local mirror lands)")
    args = ap.parse_args(argv)

    reg = _default_registry()
    if args.cmd == "status":
        for s in reg.status():
            v = f" (v{s.version})" if s.version else ""
            print(f"{s.name:16s} {s.health.value:12s}{v} {s.note or ''}")
        return 0
    if args.cmd == "refresh":
        print("refresh: not yet implemented in scaffold")
        return 0
    if args.cmd == "report":
        from .ingest.beta_matrix import read_beta_matrix
        from biocore.report.render import render_html, to_pdf
        sample = read_beta_matrix(args.sample)
        # NOTE: the live-per-marker annotation path issues one API call per
        # marker, which does not scale to a full array (935K probes = 935K
        # requests). Until the mirror-backed local lookup lands, cap the number
        # of markers annotated live so `report` on a full array cannot hang.
        # Clocks below run on the FULL sample regardless (local computation).
        markers = sample.markers[:args.max_markers]
        rep = reg.annotate(markers)
        findings = rep.all_findings()
        # epigenetic clocks: local computation, keyed by base probe id
        from . import clocks
        from .normalize import base_probe
        from biocore.providers.base import Finding, Tier, Category
        base_betas = {base_probe(k): v for k, v in sample.betas.items()}
        for cr in clocks.run_all(base_betas):
            if cr.age is None:
                continue
            findings.append(Finding(
                marker=cr.clock, source="epigenetic_clock",
                description=f"Estimated DNAm age ({cr.clock}): {cr.age:.1f} years",
                tier=(Tier.SPECULATIVE if cr.low_coverage else Tier.ROBUST),
                categories=[Category.AGING],
                detail={"coverage": cr.note}))
        html_str = render_html(findings, rep.provider_status)
        with open(args.out, "w") as fh:
            fh.write(html_str)
        print(f"{len(sample.markers)} markers -> {len(rep.all_findings())} findings -> {args.out}")
        if args.pdf:
            pdf_path = args.out.rsplit(".", 1)[0] + ".pdf"
            to_pdf(html_str, pdf_path)
            print(f"PDF: {pdf_path}")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
