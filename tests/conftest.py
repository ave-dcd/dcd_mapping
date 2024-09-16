"""Provide shared testing utilities.


Notes on test cases:
-------------------

* urn:mavedb:00000041-a-1: SRC, protein-coding, dna, uniprot ref
* urn:mavedb:00000018-a-1: HBB promoter, regulatory, DNA
* urn:mavedb:00000001-a-4: UBE2I, protein-coding, dna, uniprot ref
* urn:mavedb:00000113-a-2: APP, protein-coding, protein sequence, uniprot ref. Not in original notebooks.
* urn:mavedb:00000098-a-1: SCN5A, protein-coding, protein sequence, uniprot ref with offset
* urn:mavedb:00000061-h-1: RAF, protein coding, DNA, uniprot ref with offset
* urn:mavedb:00000068-a-1: TP53, protein-coding, DNA
"""

import json
import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from dcd_mapping.schemas import (
    AlignmentResult,
    MappedScore,
    ScoresetMetadata,
    TxSelectResult,
)

FIXTURE_DATA_DIR = Path(__file__).parents[0].resolve() / "fixtures"


def pytest_sessionstart(session) -> None:  # noqa: ARG001
    """Initialize testing environment."""
    os.environ["MAVEDB_STORAGE_DIR"] = str(FIXTURE_DATA_DIR.absolute())


@pytest.fixture(scope="session")
def fixture_data_dir():
    """Provide test data directory."""
    return FIXTURE_DATA_DIR


@pytest.fixture(scope="module")
def scoreset_metadata_fixture(fixture_data_dir: Path):
    """Provide scoreset metadata fixtures."""
    fixture_file = fixture_data_dir / "scoreset_metadata.json"
    with fixture_file.open() as f:
        data = json.load(f)
    results = {}
    for d in data["scoreset_metadata"]:
        formatted_data = ScoresetMetadata(**d)
        results[formatted_data.urn] = formatted_data
    return results


@pytest.fixture(scope="session")
def align_result_fixture(fixture_data_dir: Path):
    """Provide fixtures for alignment results."""
    fixture_file = fixture_data_dir / "align_result.json"
    with fixture_file.open() as f:
        data = json.load(f)
    results = {}
    for urn, result in data.items():
        formatted_result = AlignmentResult(**result)
        results[urn] = formatted_result
    return results


@pytest.fixture(scope="session")
def transcript_results_fixture(fixture_data_dir: Path):
    """Provide fixtures for transcript selection results."""
    fixture_file = fixture_data_dir / "transcript_result.json"
    with fixture_file.open() as f:
        data = json.load(f)
    results = {}
    for urn, result in data.items():
        formatted_result = TxSelectResult(**result)
        results[urn] = formatted_result
    return results


@pytest.fixture(scope="session")
def mapped_scores_fixture(fixture_data_dir: Path):
    fixture_file = fixture_data_dir / "mapped_scores.json"
    with fixture_file.open() as f:
        data = json.load(f)
    results = {}
    for urn, scores in data.items():
        results[urn] = [MappedScore(**score) for score in scores]
    return results


@pytest.fixture()
def mock_seqrepo_access(mocker: MagicMock):
    """Mock SeqRepo instance.

    To add new test cases, throw some print/breakpoint statements into the original
    methods and capture results there.
    """

    def _get_sequence(
        identifier: str, start: str | None = None, end: str | None = None
    ) -> str:
        calls = {
            # 41-a-1
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 14, 15): "V",
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 103, 104): "I",
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 149, 150): "Y",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 14, 15): "R",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 103, 104): "S",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 149, 150): "E",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 283, 284): "V",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 372, 373): "I",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 418, 419): "Y",
            # 99-a-1
            ("NP_000530.1", 329, 330): "D",
            ("ga4gh:SQ.kntF3MvVMQFzU94f0QegL_ktmT5f2Dk9", 329, 330): "D",
            ("ga4gh:SQ.RtClQI6dD3uj5T-DyOr9P9_-sYR2aHJX", 989, 990): "C",
            ("ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX", 129533660, 129533661): "C",
            ("ga4gh:SQ.kntF3MvVMQFzU94f0QegL_ktmT5f2Dk9", 339, 340): "T",
            ("ga4gh:SQ.RtClQI6dD3uj5T-DyOr9P9_-sYR2aHJX", 1018, 1019): "C",
            ("ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX", 129533689, 129533690): "C",
            ("ga4gh:SQ.kntF3MvVMQFzU94f0QegL_ktmT5f2Dk9", 55, 56): "F",
            ("ga4gh:SQ.RtClQI6dD3uj5T-DyOr9P9_-sYR2aHJX", 166, 167): "T",
            ("ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX", 129528899, 129528900): "T",
            ("ga4gh:SQ.kntF3MvVMQFzU94f0QegL_ktmT5f2Dk9", 3, 4): "T",
            ("ga4gh:SQ.RtClQI6dD3uj5T-DyOr9P9_-sYR2aHJX", 10, 11): "C",
            ("ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX", 129528743, 129528744): "C",
            # 103-c-1
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 336, 337): "D",
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 352, 353): "R",
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 220, 221): "M",
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 1, 2): "A",
            # 1-b-2
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 75, 76): "T",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 225, 228): "ACT",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202210743, 202210746): "AGT",
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 7, 8): "P",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 21, 24): "CCT",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 21, 23): "CC",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202220094, 202220097): "AGG",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202220095, 202220097): "GG",
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 68, 69): "Q",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 204, 207): "CAG",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 205, 207): "AG",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202210764, 202210767): "CTG",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202210764, 202210766): "CT",
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 81, 82): "M",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 245, 246): "G",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202207312, 202207313): "C",
            # 2-a-2
            ("ga4gh:SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R", 22, 23): "H",
            ("ga4gh:SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R", 10, 11): "A",
            ("ga4gh:SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R", 13, 14): "S",
            ("ga4gh:SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R", 24, 25): "D",
            ("ga4gh:SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R", 26, 27): "T",
            ("ga4gh:SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R", 25, 26): "Q",
            ("ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7", 179, 180): "A",
            ("ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7", 182, 183): "S",
            ("ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7", 191, 192): "H",
            ("ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7", 193, 194): "D",
            ("ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7", 194, 195): "Q",
            ("ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7", 195, 196): "T",
            ("NP_001123617.1", 191, 192): "H",
            ("NP_001123617.1", 179, 180): "A",
            ("NP_001123617.1", 194, 195): "Q",
            ("NP_001123617.1", 182, 183): "S",
            ("NP_001123617.1", 193, 194): "D",
            ("NP_001123617.1", 195, 196): "T",
        }
        return calls[(identifier, start, end)]

    def _translate_sequence_identifier(
        identifier: str, namespace: str | None = None
    ) -> str:
        calls = {
            ("refseq:NP_938033.1", "ga4gh"): [
                "ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye"
            ],
            ("refseq:NP_003343.1", "ga4gh"): [
                "ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip"
            ],
            ("refseq:NC_000002.12", "ga4gh"): [
                "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g"
            ],
            ("refseq:NP_002736.3", "ga4gh"): [
                "ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr"
            ],
            ("refseq:NP_001123617.1", "ga4gh"): [
                "ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7"
            ],
        }
        return calls[(identifier, namespace)]

    def _derive_refget_accession(ac: str) -> str:
        calls = {
            "NP_938033.1": "SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye",
            "NP_002736.3": "SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr",
            "NP_000530.1": "SQ.kntF3MvVMQFzU94f0QegL_ktmT5f2Dk9",
            "NP_003343.1": "SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip",
            "NC_000002.12": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
            "NC_000003.12": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "NP_001123617.1": "SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7",
        }
        return calls[ac]

    def _translate_identifier(
        ac: str, target_namespaces: str | list[str] | None = None
    ) -> tuple[list[str], str | None]:
        calls = {
            ("GRCh38:chr3", "refseq"): (["refseq:NC_000003.12"], None),
            ("GRCh37:chr3", "refseq"): (["refseq:NC_000003.11"], None),
            ("GRCh38:chr2", "refseq"): (["refseq:NC_000002.12"], None),
            ("GRCh37:chr2", "refseq"): (["refseq:NC_000002.11"], None),
        }
        return calls[(ac, target_namespaces)]

    mock_seqrepo_access = mocker.MagicMock()
    mock_seqrepo_access.get_sequence.side_effect = _get_sequence
    mock_seqrepo_access.translate_sequence_identifier.side_effect = (
        _translate_sequence_identifier
    )
    mock_seqrepo_access.derive_refget_accession.side_effect = _derive_refget_accession
    mock_seqrepo_access.translate_identifier.side_effect = _translate_identifier
    mocker.patch("dcd_mapping.vrs_map.get_seqrepo", return_value=mock_seqrepo_access)
    mocker.patch("dcd_mapping.lookup.get_seqrepo", return_value=mock_seqrepo_access)
    mocker.patch("dcd_mapping.annotate.get_seqrepo", return_value=mock_seqrepo_access)
    return mock_seqrepo_access
