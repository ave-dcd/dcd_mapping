"""Provide shared testing utilities."""
import json
from pathlib import Path

import pytest

from dcd_mapping.schemas import AlignmentResult, ScoresetMetadata, TxSelectResult


@pytest.fixture(scope="module")
def scoreset_metadata_fixture():
    """Provide scoreset metadata fixtures."""
    fixture_file = (
        Path(__file__).parents[0].resolve() / "fixtures" / "scoreset_metadata.json"
    )
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for d in data["scoreset_metadata"]:
        formatted_data = ScoresetMetadata(**d)
        results[formatted_data.urn] = formatted_data
    return results


@pytest.fixture(scope="session")
def align_result_fixture():
    """Provide fixtures for alignment results."""
    fixture_file = (
        Path(__file__).parents[0].resolve() / "fixtures" / "align_result.json"
    )
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for urn, result in data.items():
        formatted_result = AlignmentResult(**result)
        results[urn] = formatted_result
    return results


@pytest.fixture(scope="session")
def transcript_results_fixture():
    """Provide fixtures for transcript selection results."""
    fixture_file = (
        Path(__file__).parents[0].resolve() / "fixtures" / "transcript_result.json"
    )
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for urn, result in data.items():
        formatted_result = TxSelectResult(**result)
        results[urn] = formatted_result
    return results
