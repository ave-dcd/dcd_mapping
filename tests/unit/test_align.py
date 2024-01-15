"""Test ``align`` module.

Todo:
----
* Mock the BLAT call/result file

"""
import pytest
from cool_seq_tool.schemas import Strand

from dcd_mapping.align import align


def test_align_src_catalytic_domain(scoreset_metadata_fixture):
    """Test ``align()`` method on urn:mavedb:00000041-a-1"""
    scoreset_metadata = scoreset_metadata_fixture["urn:mavedb:00000041-a-1"]
    align_result = align(scoreset_metadata)
    assert align_result
    assert align_result.chrom == "chr20"
    assert align_result.strand == Strand.POSITIVE
    assert align_result.coverage == pytest.approx(100.0)
    assert align_result.ident_pct == pytest.approx(99.86)
    assert align_result.query_range.start == 0
    assert align_result.query_range.end == 750
    query_subranges = [
        [0, 52],
        [52, 232],
        [232, 309],
        [309, 463],
        [463, 595],
        [595, 750],
    ]
    for actual, expected in zip(align_result.query_subranges, query_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.query_subranges) == len(query_subranges)
    assert align_result.hit_range.start == 37397802
    assert align_result.hit_range.end == 37403325
    hit_subranges = [
        [37397802, 37397854],
        [37400114, 37400294],
        [37401601, 37401678],
        [37402434, 37402588],
        [37402748, 37402880],
        [37403170, 37403325],
    ]
    for actual, expected in zip(align_result.hit_subranges, hit_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.hit_subranges) == len(hit_subranges)


def test_align_hbb(scoreset_metadata_fixture):
    """Test ``align()`` method on urn:mavedb:00000018-a-1"""
    scoreset_metadata = scoreset_metadata_fixture["urn:mavedb:00000018-a-1"]
    align_result = align(scoreset_metadata)
    assert align_result
    assert align_result.chrom == "chr11"
    assert align_result.strand == Strand.POSITIVE
    assert align_result.coverage == pytest.approx(100.0)
    assert align_result.ident_pct == 0  # TODO
    assert align_result.query_range.start == 0
    assert align_result.query_range.end == 187
    query_subranges = [[0, 187]]
    for actual, expected in zip(align_result.query_subranges, query_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.query_subranges) == len(query_subranges)
    assert align_result.hit_range.start == 5227021
    assert align_result.hit_range.end == 5227208
    hit_subranges = [[5227021, 5227208]]
    for actual, expected in zip(align_result.hit_subranges, hit_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.hit_subranges) == len(hit_subranges)


def test_align_ube2i(scoreset_metadata_fixture):
    """Test ``align()`` on urn:mavedb:00000001-a-4"""
    scoreset_metadata = scoreset_metadata_fixture["urn:mavedb:00000001-a-4"]
    align_result = align(scoreset_metadata)
    assert align_result
    assert align_result.chrom == "chr16"
    assert align_result.strand == Strand.POSITIVE
    assert align_result.coverage == pytest.approx(100.0)
    assert align_result.ident_pct == pytest.approx(100.0)
    assert align_result.query_range.start == 0
    assert align_result.query_range.end == 477
    query_subranges = [
        [0, 66],
        [66, 150],
        [150, 223],
        [223, 333],
        [333, 413],
        [413, 477],
    ]
    for actual, expected in zip(align_result.query_subranges, query_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.query_subranges) == len(query_subranges)
    assert align_result.hit_range.start == 5227021
    assert align_result.hit_range.end == 5227208
    hit_subranges = [
        [1314030, 1314096],
        [1314292, 1314376],
        [1315653, 1315726],
        [1320173, 1320283],
        [1320437, 1320517],
        [1324729, 1324793],
    ]
    for actual, expected in zip(align_result.hit_subranges, hit_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.hit_subranges) == len(hit_subranges)


def test_align_scn5a(scoreset_metadata_fixture):
    """Test ``align()`` method on urn:mavedb:00000098-a-1"""
    scoreset_metadata = scoreset_metadata_fixture["urn:mavedb:00000098-a-1"]
    align_result = align(scoreset_metadata)
    assert align_result
    assert align_result.chrom == "chr3"
    assert align_result.strand == Strand.NEGATIVE
    assert align_result.coverage == pytest.approx(100.0)
    assert align_result.ident_pct == pytest.approx(100.0)
    assert align_result.query_range.start == 0
    assert align_result.query_range.end == 12
    query_subranges = [[0, 12]]
    for actual, expected in zip(align_result.query_subranges, query_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.query_subranges) == len(query_subranges)
    assert align_result.hit_range.start == 38551475
    assert align_result.hit_range.end == 38551511
    hit_subranges = [[38551475, 38551511]]
    for actual, expected in zip(align_result.hit_subranges, hit_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.hit_subranges) == len(hit_subranges)


def test_align_raf(scoreset_metadata_fixture):
    """Test ``align()`` method on urn:mavedb:00000061-h-1."""
    scoreset_metadata = scoreset_metadata_fixture["urn:mavedb:00000061-h-1"]
    align_result = align(scoreset_metadata)
    assert align_result
    assert align_result.chrom == "chr3"
    assert align_result.strand == Strand.NEGATIVE
    assert align_result.coverage == pytest.approx(100.0)
    assert align_result.ident_pct == pytest.approx(100.0)
    assert align_result.query_range.start == 0
    assert align_result.query_range.end == 117
    query_subranges = [[54, 117], [0, 54]]
    for actual, expected in zip(align_result.query_subranges, query_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.query_subranges) == len(query_subranges)
    assert align_result.hit_range.start == 12618514
    assert align_result.hit_range.end == 12612062
    hit_subranges = [[12611999, 12612062], [12618514, 12618568]]
    for actual, expected in zip(align_result.hit_subranges, hit_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert len(align_result.hit_subranges) == len(hit_subranges)
