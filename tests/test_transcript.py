"""Test ``transcripts`` module.

Todo:
----
* Get a test case where there are no common transcripts and you have to use the UniProt ref.


"""

import re
from collections.abc import Coroutine
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pytest

from dcd_mapping.mavedb_data import _load_scoreset_records, get_scoreset_records
from dcd_mapping.schemas import AlignmentResult, ScoresetMetadata, TxSelectResult
from dcd_mapping.transcripts import select_transcript


@pytest.fixture()
def mock_cst(mocker: MagicMock, mock_seqrepo_access):
    """Mock CoolSeqTool instance."""

    async def _execute_query(query: str) -> Coroutine[Any, Any, Any]:
        query = query.strip().replace("\n", "")
        query = re.sub(r"\s+", " ", query)
        calls = {
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202207252 BETWEEN alt_start_i AND alt_end_i OR 202207321 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005782.1"},
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_001005782.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371393.1"},
                {"tx_ac": "NM_001371394.1"},
                {"tx_ac": "NM_003352.8"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202210734 BETWEEN alt_start_i AND alt_end_i OR 202210806 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_001005782.1"},
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005782.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371394.1"},
                {"tx_ac": "NM_003352.8"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202214356 BETWEEN alt_start_i AND alt_end_i OR 202214434 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005782.1"},
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_001005782.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371393.1"},
                {"tx_ac": "NM_001371394.1"},
                {"tx_ac": "NM_003352.8"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SUMO1' AND (202220031 BETWEEN alt_start_i AND alt_end_i OR 202220109 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000002.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_001005781.1"},
                {"tx_ac": "NM_003352.8"},
                {"tx_ac": "NM_003352.4"},
                {"tx_ac": "NM_001005781.2"},
                {"tx_ac": "NM_001371392.1"},
                {"tx_ac": "NM_001371393.1"},
                {"tx_ac": "NM_001371394.1"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37397802 BETWEEN alt_start_i AND alt_end_i OR 37397854 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37400114 BETWEEN alt_start_i AND alt_end_i OR 37400294 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37401601 BETWEEN alt_start_i AND alt_end_i OR 37401678 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%'; ": [
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37402434 BETWEEN alt_start_i AND alt_end_i OR 37402588 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37402748 BETWEEN alt_start_i AND alt_end_i OR 37402880 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SRC' AND (37403170 BETWEEN alt_start_i AND alt_end_i OR 37403325 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000020.11' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_198291.2"},
                {"tx_ac": "NM_005417.5"},
                {"tx_ac": "NM_005417.4"},
                {"tx_ac": "NM_198291.3"},
            ],
            "SELECT tx_ac FROM uta_20210129b.tx_exon_aln_v WHERE hgnc = 'SCN5A' AND (38551475 BETWEEN alt_start_i AND alt_end_i OR 38551511 BETWEEN alt_start_i AND alt_end_i) AND alt_ac = 'NC_000003.12' AND tx_ac NOT LIKE 'NR_%';": [
                {"tx_ac": "NM_198056.2"},
                {"tx_ac": "NM_001160161.1"},
                {"tx_ac": "NM_001354701.1"},
                {"tx_ac": "NM_001099405.1"},
                {"tx_ac": "NM_000335.4"},
                {"tx_ac": "NM_001160160.1"},
                {"tx_ac": "NM_001099404.1"},
                {"tx_ac": "NM_001160160.2"},
                {"tx_ac": "NM_001160161.2"},
                {"tx_ac": "NM_001354701.2"},
                {"tx_ac": "NM_001099404.2"},
                {"tx_ac": "NM_001099405.2"},
                {"tx_ac": "NM_000335.5"},
                {"tx_ac": "NM_198056.3"},
            ],
        }
        return calls[query]

    mock_uta_instance = mocker.MagicMock()
    mock_uta_instance.schema = "20210129b"
    mock_uta_instance.execute_query.side_effect = _execute_query

    mock_cst_instance = mocker.MagicMock()
    mock_cst_instance.uta = mock_uta_instance
    mock_cst_instance.seqrepo_access = mock_seqrepo_access
    mocker.patch(
        "dcd_mapping.lookup.CoolSeqToolBuilder.__new__", return_value=mock_cst_instance
    )

    return mock_cst_instance


def check_transcript_results_equality(actual: TxSelectResult, expected: TxSelectResult):
    """Check equality of transcript selection result vs fixture"""
    assert actual.np == expected.np
    assert actual.start == expected.start
    assert actual.is_full_match is expected.is_full_match
    assert actual.nm == expected.nm
    assert actual.transcript_mode == expected.transcript_mode


@pytest.mark.asyncio(scope="module")
async def test_1_b_2(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
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
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
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
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
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
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
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
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
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
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
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
