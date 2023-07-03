# TODO: create a separate directory for test
import pytest
from blat_alignment import mave_to_blat, check_non_human
from metadata_process import metadata_obtain


class TestBLATHuman:
    def test_for_human_urn(self):
        with open("urn:mavedb:00000041-b-1") as scoreset:
            mave_dat = metadata_obtain(scoreset)

        mave_blat = mave_to_blat(mave_dat)
        organism = check_non_human(mave_blat)
        assert organism == "human"

    def test_for_non_human_urn(self):
        with open("urn:mavedb:00000004-a-1") as scoreset:
            mave_dat = metadata_obtain(scoreset)

        mave_blat = mave_to_blat(mave_dat)
        organism = check_non_human(mave_blat)
        assert organism == "Non human"

    def test_for_human_urn_2(self):
        with open("urn:mavedb:00000083-e-1") as scoreset:
            mave_dat = metadata_obtain(scoreset)

        mave_blat = mave_to_blat(mave_dat)
        organism = check_non_human(mave_blat)
        assert organism == "human"

    def test_for_non_human_urn_2(self):
        with open("urn:mavedb:00000010-a-1") as scoreset:
            mave_dat = metadata_obtain(scoreset)

        mave_blat = mave_to_blat(mave_dat)
        organism = check_non_human(mave_blat)
        assert organism == "Non human"
