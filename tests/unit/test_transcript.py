"""Test ``transcripts`` module."""
from typing import Dict

import pytest

from dcd_mapping.resources import get_scoreset_records
from dcd_mapping.schemas import AlignmentResult, ScoresetMetadata, TxSelectResult
from dcd_mapping.transcripts import select_transcript


@pytest.mark.asyncio(scope="module")
async def test_tx_src(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    """Test transcript selection for urn:mavedb:00000041-a-1"""
    urn = "urn:mavedb:00000041-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)  # TODO real fixture
    alignment_result = align_result_fixture[urn]
    expected = transcript_results_fixture[urn]

    actual = await select_transcript(metadata, records, alignment_result)

    assert actual
    assert actual.np == expected.np
    assert actual.start == expected.start
    assert actual.is_full_match is expected.is_full_match
    assert actual.nm == expected.nm
    assert actual.transcript_mode == expected.transcript_mode


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
    assert actual.np == expected.np
    assert actual.start == expected.start
    assert actual.is_full_match is expected.is_full_match
    assert actual.nm == expected.nm
    assert actual.transcript_mode == expected.transcript_mode


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
    assert actual.np == expected.np
    assert actual.start == expected.start
    assert actual.is_full_match is expected.is_full_match
    assert actual.nm == expected.nm
    assert actual.transcript_mode == expected.transcript_mode
