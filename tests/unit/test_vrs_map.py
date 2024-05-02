"""Test ``vrs_map.py``

Todo:
----
* Sample a second mapping row for each test case (also, order shouldn't matter)
* Add a few more test cases


"""
from pathlib import Path
from typing import Dict

import pytest

from dcd_mapping.mavedb_data import _load_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    TxSelectResult,
)
from dcd_mapping.vrs_map import vrs_map


@pytest.fixture()
def mock_seqrepo(mocker):
    mock_seqrepo_instance = mocker.MagicMock()

    def _get_sequence(identifier, start=None, end=None):
        calls = {
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 14, 15): "V",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 283, 284): "V",
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 149, 150): "Y",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 418, 419): "Y",
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 103, 104): "I",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 372, 373): "I",
        }
        return calls[(identifier, start, end)]

    mock_seqrepo_access = mocker.MagicMock()
    mock_seqrepo_access.get_sequence.side_effect = _get_sequence
    mock_seqrepo_access.sr = mock_seqrepo_instance
    mocker.patch("dcd_mapping.vrs_map.get_seqrepo", return_value=mock_seqrepo_access)
    return mock_seqrepo_access


def test_41_a_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
    mock_seqrepo,
):
    urn = "urn:mavedb:00000041-a-1"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_result = transcript_results_fixture[urn]

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None
    pre_mapped = mappings[0].pre_mapped_variants
    assert pre_mapped["id"] == "ga4gh:VA.NJgaCF0JPFERdw9Y7fW4bXkaP1tSa5fv"
    assert pre_mapped["location"]["interval"]["start"]["value"] == 14
    assert pre_mapped["location"]["interval"]["end"]["value"] == 15
    assert pre_mapped["state"]["sequence"] == "T"
    assert (
        pre_mapped["location"]["sequence_id"]
        == "ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB"
    )
    post_mapped = mappings[0].post_mapped_variants
    assert post_mapped["id"] == "ga4gh:VA.csCB31gWoiiD38TlR35dZnyAI156YWgW"
    assert post_mapped["location"]["interval"]["start"]["value"] == 283
    assert post_mapped["location"]["interval"]["end"]["value"] == 284
    assert post_mapped["state"]["sequence"] == "T"
    assert (
        post_mapped["location"]["sequence_id"]
        == "ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye"
    )

    mock_seqrepo.sr.store.assert_any_call(
        "LRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGETGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKEPEERPTFEYLQAFL",
        nsaliases=[
            {"namespace": "ga4gh", "alias": "SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB"}
        ],
    )
    mock_seqrepo.sr.store.assert_any_call(
        "CTGCGGCTGGAGGTCAAGCTGGGCCAGGGCTGCTTTGGCGAGGTGTGGATGGGGACCTGGAACGGTACCACCAGGGTGGCCATCAAAACCCTGAAGCCTGGCACGATGTCTCCAGAGGCCTTCCTGCAGGAGGCCCAGGTCATGAAGAAGCTGAGGCATGAGAAGCTGGTGCAGTTGTATGCTGTGGTTTCAGAGGAGCCCATTTACATCGTCACGGAGTACATGAGCAAGGGGAGTTTGCTGGACTTTCTCAAGGGGGAGACAGGCAAGTACCTGCGGCTGCCTCAGCTGGTGGACATGGCTGCTCAGATCGCCTCAGGCATGGCGTACGTGGAGCGGATGAACTACGTCCACCGGGACCTTCGTGCAGCCAACATCCTGGTGGGAGAGAACCTGGTGTGCAAAGTGGCCGACTTTGGGCTGGCTCGGCTCATTGAAGACAATGAGTACACGGCGCGGCAAGGTGCCAAATTCCCCATCAAGTGGACGGCTCCAGAAGCTGCCCTCTATGGCCGCTTCACCATCAAGTCGGACGTGTGGTCCTTCGGGATCCTGCTGACTGAGCTCACCACAAAGGGACGGGTGCCCTACCCTGGGATGGTGAACCGCGAGGTGCTGGACCAGGTGGAGCGGGGCTACCGGATGCCCTGCCCGCCGGAGTGTCCCGAGTCCCTGCACGACCTCATGTGCCAGTGCTGGCGGAAGGAGCCTGAGGAGCGGCCCACCTTCGAGTACCTGCAGGCCTTCCTG",
        nsaliases=[
            {"namespace": "ga4gh", "alias": "SQ.SBlyhKAQ5u14huYar7I7UmXZk4eRPexx"}
        ],
    )
