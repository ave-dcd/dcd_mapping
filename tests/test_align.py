"""Test ``align`` module.

Todo:
----
* Mock the BLAT call/result file


"""

import pytest

from dcd_mapping.align import align
from dcd_mapping.schemas import AlignmentResult


def check_alignment_result_equality(actual: AlignmentResult, expected: AlignmentResult):
    """Check correctness of alignment result against fixture"""
    assert actual.chrom == expected.chrom
    assert actual.strand == expected.strand
    assert actual.coverage == pytest.approx(expected.coverage)
    assert actual.ident_pct == pytest.approx(expected.ident_pct)
    assert actual.query_range.start == expected.query_range.start
    assert actual.query_range.end == expected.query_range.end
    for a, e in zip(actual.query_subranges, expected.query_subranges, strict=False):
        assert a.start == e.start
        assert a.end == e.end
    assert len(actual.query_subranges) == len(expected.query_subranges)
    assert actual.hit_range.start == expected.hit_range.start
    assert actual.hit_range.end == expected.hit_range.end
    for a, e in zip(actual.hit_subranges, expected.hit_subranges, strict=False):
        assert a.start == e.start
        assert a.end == e.end
    assert len(actual.hit_subranges) == len(expected.hit_subranges)


def test_align_src_catalytic_domain(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000041-a-1"""
    urn = "urn:mavedb:00000041-a-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = align_result_fixture[urn]
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_hbb(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000018-a-1"""
    urn = "urn:mavedb:00000018-a-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = align_result_fixture[urn]
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_ube2i(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` on urn:mavedb:00000001-a-4"""
    urn = "urn:mavedb:00000001-a-4"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = align_result_fixture[urn]
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_scn5a(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000098-a-1"""
    urn = "urn:mavedb:00000098-a-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = align_result_fixture[urn]
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_raf(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000061-h-1."""
    urn = "urn:mavedb:00000061-h-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = align_result_fixture[urn]
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_tp53(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000068-a-1"""
    urn = "urn:mavedb:00000068-a-1"
    metadata = scoreset_metadata_fixture[urn]
    align_result = align(metadata)
    expected = align_result_fixture[urn]
    assert align_result
    check_alignment_result_equality(align_result, expected)
