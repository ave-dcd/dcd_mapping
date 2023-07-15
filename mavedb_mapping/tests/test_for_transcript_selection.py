import pytest
import pickle
from metadata_process import metadata_obtain
from blat_alignment import mave_to_blat
from transcript_selection import main

file = open("mappings.pickle", "rb")
mappings = pickle.load(file)
file.close()


def test_for_refseq_prot():
    with open("tests/data/urn:mavedb:00000041-b-1") as scoreset:
        mave_dat = metadata_obtain(scoreset)
    mave = mave_to_blat(mave_dat)
    tr = main(mave, mave_dat)
    computed = tr["RefSeq_prot"]
    assert computed == mappings["urn:mavedb:00000041-a-1"][0]


def test_for_refseq_nuc():
    with open("tests/data/urn:mavedb:00000060-a-1") as scoreset:
        mave_dat = metadata_obtain(scoreset)
    mave = mave_to_blat(mave_dat)
    tr = main(mave, mave_dat)
    computed = tr["RefSeq_nuc"]
    assert computed == mappings["urn:mavedb:00000060-a-1"][4]
