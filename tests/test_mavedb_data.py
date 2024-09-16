"""Test `dcd_mapping.resources` module."""

import json
import shutil
from pathlib import Path

import pytest
import requests_mock

from dcd_mapping.mavedb_data import get_scoreset_metadata, get_scoreset_records


@pytest.fixture()
def resources_data_dir():
    """Temporarily store data resources"""
    path = Path(__file__).parent / "tmp"
    if path.exists():  # make sure it's empty
        shutil.rmtree(str(path.absolute()))
    else:
        path.mkdir()
    yield path
    shutil.rmtree(str(path.absolute()))  # clean up afterward


@pytest.fixture()
def scoreset_metadata_response(fixture_data_dir: Path):
    """Provide response that client receives from MaveDB API"""
    with (fixture_data_dir / "scoreset_metadata_response.json").open() as f:
        return json.load(f)


def test_get_scoreset_metadata(
    resources_data_dir: Path, scoreset_metadata_response: dict
):
    urn = "urn:mavedb:00000093-a-1"
    with requests_mock.Mocker() as m:
        m.get(
            f"https://api.mavedb.org/api/v1/score-sets/{urn}",
            json=scoreset_metadata_response[urn],
        )
        scoreset_metadata = get_scoreset_metadata(
            urn, dcd_mapping_dir=resources_data_dir
        )
        assert scoreset_metadata.urn == urn
        assert scoreset_metadata.target_uniprot_ref
        assert scoreset_metadata.target_uniprot_ref.id == "uniprot:P38398"
        assert scoreset_metadata.target_uniprot_ref.offset == 0


def test_get_scoreset_records(resources_data_dir: Path, fixture_data_dir: Path):
    urn = "urn:mavedb:00000093-a-1"
    with (fixture_data_dir / f"{urn}_scores.csv").open() as f:
        scores_csv_text = f.read()
    with requests_mock.Mocker() as m:
        m.get(
            f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores",
            text=scores_csv_text,
        )
        scoreset_records = get_scoreset_records(urn, dcd_mapping_dir=resources_data_dir)
        assert len(scoreset_records) == 853
