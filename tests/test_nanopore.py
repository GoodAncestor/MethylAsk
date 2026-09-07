"""Synthetic native-ONT mapping regressions; no human sample or network access."""
import gzip
import pytest

from methylask.ingest.nanopore import read_bedmethyl, read_modbam

MAP = {"cg00000001": ("chr1", 10), "cg00000002": ("chr1", 20),
       "cg00000003": ("chr1", 30)}


def row(pos=10, mod="m", strand="+", nmod=2, ncan=2, other=1, chrom="1"):
    cov = nmod + ncan + other
    return "\t".join(map(str, [chrom, pos, pos+1, mod, cov, strand, pos, pos+1,
                              "0,0,0", cov, 100*nmod/cov if cov else 0,
                              nmod, ncan, other, 0, 0, 0, 0])) + "\n"


def bed(tmp_path, text):
    path = tmp_path / "calls.bed.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return str(path)


def read(path, **kwargs):
    return read_bedmethyl(path, reference_build="GRCh38", probe_map=MAP, **kwargs)


def test_counts_strands_modifications_and_missingness(tmp_path):
    path = bed(tmp_path, row(nmod=2, ncan=0, other=1) + row(pos=11, strand="-", nmod=1, ncan=1, other=0)
               + row(mod="h", nmod=1, ncan=0, other=2) + row(pos=20, nmod=0, ncan=1, other=0)
               + row(pos=99))
    sample = read(path)
    assert sample.betas == {"cg00000001": 3/5}
    assert sample.coverage == {"cg00000001": 5, "cg00000002": 1}
    assert sample.stats["missing_probes"] == 1
    assert sample.stats["low_coverage_probes"] == 1
    assert sample.stats["mapped_probes"] == 1
    assert sample.stats["rows"] == 5
    assert sample.provenance["mod_code"] == "m"
    assert sample.provenance["array_clocks_allowed_by_default"] is False
    assert sample.warnings and "5hmC" in sample.warnings[0]


def test_hydroxymethylation_is_selectable_and_distinct(tmp_path):
    path = bed(tmp_path, row(mod="m") + row(mod="h", nmod=1, ncan=2, other=2))
    assert read(path).betas["cg00000001"] == .4
    assert read(path, mod_code="h").betas["cg00000001"] == .2


def test_plain_and_annotated_names_spaces_and_combined_strands(tmp_path):
    path = bed(tmp_path, row(mod="m,CG,0", strand=".").replace("\t", " "))
    with pytest.raises(ValueError, match="combined_strands"):
        read(path)
    assert read(path, combined_strands=True).betas == {"cg00000001": .4}
    with pytest.raises(ValueError, match="exclusively"):
        read(bed(tmp_path, row()), combined_strands=True)


def test_duplicate_motif_rows_cannot_inflate_coverage(tmp_path):
    with pytest.raises(ValueError, match="Duplicate"):
        read(bed(tmp_path, row() + row(mod="m,CG,0")))


@pytest.mark.parametrize("text", ["1\t10\t11\tm\n", row().replace("\t2\t2\t1\t", "\t9\t2\t1\t"),
                                  row().replace("1\t10\t11", "1\t-1\t0")])
def test_malformed_input_is_not_silent_empty_report(tmp_path, text):
    with pytest.raises(ValueError, match="bedMethyl row"):
        read(bed(tmp_path, text))


def test_empty_and_uncovered_are_explicit(tmp_path):
    sample = read(bed(tmp_path, row(nmod=0, ncan=0, other=0)))
    assert sample.betas == {}
    assert sample.coverage == {"cg00000001": 0}
    assert sample.stats["status"] == "no_qualifying_markers"
    with pytest.raises(ValueError, match="GRCh38"):
        read_bedmethyl("unused", reference_build="hg19", probe_map=MAP)
    with pytest.raises(ValueError, match="No unambiguous"):
        read_bedmethyl("unused", reference_build="GRCh38", probe_map={})


def test_manifest_replica_ids_deduplicated_or_excluded(tmp_path):
    mapping = {"cg00000001_A": ("chr1", 10), "cg00000001_B": ("1", 10),
               "cg00000002_A": ("chr1", 20), "cg00000002_B": ("chr1", 30)}
    sample = read_bedmethyl(bed(tmp_path, row()), reference_build="hg38", probe_map=mapping)
    assert sample.betas == {"cg00000001": .4}
    assert sample.stats["total_probes"] == 1
    assert sample.stats["ambiguous_manifest_probes"] == 1


def bam(tmp_path, records, *, sample_ids=("synthetic",)):
    pysam = pytest.importorskip("pysam")
    ref = "A" * 10 + "CG" + "A" * 28
    fa = str(tmp_path / "ref.fa")
    with open(fa, "w") as fh:
        fh.write(">1\n" + ref + "\n")
    pysam.faidx(fa)
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": len(ref)}],
              "RG": [{"ID": str(i), "SM": s} for i, s in enumerate(sample_ids)]}
    path = str(tmp_path / "calls.bam")
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i, (flag, mm, ml) in enumerate(records):
            r = pysam.AlignedSegment()
            r.query_name = f"synthetic{i}"
            r.query_sequence = "CG"
            r.flag = flag
            r.reference_id = 0
            r.reference_start = 10
            r.cigar = [(0, 2)]
            r.mapping_quality = 60
            r.set_tag("RG", "0")
            if mm is not None:
                r.set_tag("MM", mm)
            if ml is not None:
                r.set_tag("ML", ml)
            out.write(r)
    return path, fa


def test_modbam_reverse_other_mod_ambiguous_implicit_and_filtering(tmp_path):
    records = [(0, "C+mh?,0;", [240, 1]), (16, "C+mh?,0;", [240, 1]),
               (0, "C+mh?,0;", [1, 240]),  # h belongs only in denominator
               (0, "C+mh?,0;", [1, 1]),    # canonical
               (0, "C+mh.;", None),        # implicit canonical
               (0, "C+mh?,0;", [120, 120]), # ambiguous -> not canonical
               (0, "C+mh?;", None),        # uncalled -> not canonical
               (1024, "C+mh?,0;", [240, 1]),
               (256, "C+mh?,0;", [240, 1]),
               (2048, "C+mh?,0;", [240, 1]), (0, None, None)]
    path, fa = bam(tmp_path, records)
    sample = read_modbam(path, reference_build="GRCh38", reference_fasta=fa, probe_map=MAP)
    assert sample.betas == {"cg00000001": 2/5}
    assert sample.coverage["cg00000001"] == 5
    assert sample.stats["ambiguous_calls"] == 1
    assert sample.stats["unknown_calls"] == 1
    assert sample.stats["reads_filtered"] == 3
    assert sample.provenance["sample_ids"] == ["synthetic"]


def test_modbam_default_mm_flag_and_reverse_implicit(tmp_path):
    path, fa = bam(tmp_path, [(0, "C+m;", None), (16, "C+m.;", None)])
    sample = read_modbam(path, reference_build="hg38", reference_fasta=fa,
                         probe_map=MAP, min_coverage=2)
    assert sample.betas == {"cg00000001": 0.0}


def test_modbam_multiple_people_fail(tmp_path):
    path, fa = bam(tmp_path, [(0, "C+m;", None)], sample_ids=("one", "two"))
    with pytest.raises(ValueError, match="multiple RG SM"):
        read_modbam(path, reference_build="hg38", reference_fasta=fa, probe_map=MAP)


def test_modbam_unknown_probability_not_canonical(tmp_path):
    path, fa = bam(tmp_path, [(0, "C+m?,0;", None)])
    sample = read_modbam(path, reference_build="hg38", reference_fasta=fa, probe_map=MAP)
    assert sample.betas == {} and sample.stats["unknown_calls"] == 1


def test_modbam_and_bedmethyl_agree_on_counts(tmp_path):
    path, fa = bam(tmp_path, [(0, "C+mh?,0;", [240, 1]), (16, "C+mh?,0;", [1, 240])])
    sample = read_modbam(path, reference_build="GRCh38", reference_fasta=fa, probe_map=MAP, min_coverage=2)
    other = read(bed(tmp_path, row(nmod=1, ncan=0, other=0)
                     + row(pos=11, strand="-", nmod=0, ncan=0, other=1)), min_coverage=2)
    assert sample.betas == other.betas == {"cg00000001": .5}
    assert sample.coverage == other.coverage
