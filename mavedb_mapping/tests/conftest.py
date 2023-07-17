import pytest
from blat_alignment import mave_to_blat
from metadata_process import metadata_obtain
from transcript_selection import main
from vrs_mapping import vrs_mapping


@pytest.fixture(
    scope="package", params=["urn:mavedb:00000041-a-1", "urn:mavedb:00000001-a-4"]
)
def full_mapping(request):
    scoreset_path = f"tests/data/{request.param}"
    with open(scoreset_path) as scoreset:
        mave_dat = metadata_obtain(scoreset)
    mave_blat = mave_to_blat(mave_dat)
    mappings_dict = main(mave_blat, mave_dat)
    vrs_mapped = vrs_mapping(mave_dat, mappings_dict, mave_blat)
    return vrs_mapped, mave_dat["urn"]
