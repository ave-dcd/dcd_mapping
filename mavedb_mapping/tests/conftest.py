import pytest
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.transcript_selection import main
from mavedb_mapping.vrs_mapping import vrs_mapping
from mavedb_mapping import data_file_path


@pytest.fixture(
    scope="package", params=["urn:mavedb:00000041-a-1", "urn:mavedb:00000001-a-4"]
)
def full_mapping(request):
    scoreset_path = f"{data_file_path}{request.param}"
    with open(scoreset_path) as scoreset:
        mave_dat = metadata_obtain(scoreset)
    mave_blat = mave_to_blat(mave_dat)
    mappings_dict = main(mave_blat, mave_dat)
    vrs_mapped = vrs_mapping(mave_dat, mappings_dict, mave_blat)
    return vrs_mapped, mave_dat["urn"]
