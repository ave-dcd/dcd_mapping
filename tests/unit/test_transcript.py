"""Test ``transcripts`` module."""
from pathlib import Path
from typing import Dict

import pytest

from dcd_mapping.mavedb_data import _load_scoreset_records, get_scoreset_records
from dcd_mapping.schemas import AlignmentResult, ScoresetMetadata, TxSelectResult
from dcd_mapping.transcripts import select_transcript


def check_transcript_results_equality(actual: TxSelectResult, expected: TxSelectResult):
    """Check equality of transcript selection result vs fixture"""
    assert actual.np == expected.np
    assert actual.start == expected.start
    assert actual.is_full_match is expected.is_full_match
    assert actual.nm == expected.nm
    assert actual.transcript_mode == expected.transcript_mode


@pytest.mark.asyncio()
async def test_1_b_2(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    urn = "urn:mavedb:00000001-b-2"
    metadata = scoreset_metadata_fixture[urn]
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    align_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, align_result)
    assert actual, "`select_transcript()` should return a transcript selection result"
    check_transcript_results_equality(actual, expected)


@pytest.mark.asyncio(scope="module")
async def test_tx_src(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000041-a-1"""
    urn = "urn:mavedb:00000041-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)


@pytest.mark.asyncio(scope="module")
async def test_tx_scn5a(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000098-a-1"""
    urn = "urn:mavedb:00000098-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)


@pytest.mark.asyncio(scope="module")
async def test_tx_hbb(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000018-a-1"""
    urn = "urn:mavedb:00000018-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual is None


@pytest.mark.asyncio(scope="module")
async def test_tx_raf(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000061-h-1"""
    urn = "urn:mavedb:00000061-h-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)


@pytest.mark.asyncio(scope="module")
async def test_tx_tp53(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000068-a-1"""
    urn = "urn:mavedb:00000068-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]
    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    check_transcript_results_equality(actual, expected)
