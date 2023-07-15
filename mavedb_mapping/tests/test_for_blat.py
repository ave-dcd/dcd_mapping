# TODO: create a separate directory for test
import pytest
from blat_alignment import mave_to_blat, check_non_human
from metadata_process import metadata_obtain


@pytest.fixture
def scoreset_organism(request):
    scoreset_path = request.param
    with open(scoreset_path) as scoreset:
        mave_dat = metadata_obtain(scoreset)
    mave_blat = mave_to_blat(mave_dat)
    return mave_blat


@pytest.mark.parametrize(
    "scoreset_organism",
    ["tests/data/urn:mavedb:00000041-b-1", "tests/data/urn:mavedb:00000083-e-1"],
    indirect=True,
)
def test_human_organism(scoreset_organism):
    organism = check_non_human(scoreset_organism)
    assert organism == "human"


@pytest.mark.parametrize(
    "scoreset_organism",
    [
        "tests/data/urn:mavedb:00000010-a-1",
        "tests/data/urn:mavedb:00000004-a-1",
    ],
    indirect=True,
)
def test_non_human_organism(scoreset_organism):
    organism = check_non_human(scoreset_organism)
    assert organism == "Non human"
