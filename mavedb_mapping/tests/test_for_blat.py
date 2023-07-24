import pytest
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.transcript_selection_helper import HelperFunctionsForBLATOutput

"""Tests that run check_non_human function after BLAT Alignment to determine if scoreset is Human"""

data_file_path = "mavedb_mapping/tests/data/"


@pytest.fixture
def blat_organism(request):
    """Fixture to return dictionary after BLAT Alignment"""
    scoreset_path = f"{data_file_path}{request.param}"
    with open(scoreset_path) as scoreset:
        mave_dat = metadata_obtain(scoreset)
    mave_blat = mave_to_blat(mave_dat)
    helper = HelperFunctionsForBLATOutput(mave_blat)
    organism = helper.is_human()
    return organism


@pytest.mark.parametrize(
    "blat_organism",
    ["urn:mavedb:00000041-b-1", "urn:mavedb:00000083-e-1"],
    indirect=True,
)
def test_human_organism(blat_organism):
    """Test to check if input scoreset organism is Human"""
    assert blat_organism == True


@pytest.mark.parametrize(
    "blat_organism",
    [
        "urn:mavedb:00000010-a-1",
        "urn:mavedb:00000004-a-1",
    ],
    indirect=True,
)
def test_non_human_organism(blat_organism):
    """Test to check if input scoreset organism is Non human"""
    assert blat_organism == False
