import pytest
from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.transcript_selection import main


@pytest.mark.parametrize(
    "obtain_transcripts", ["urn:mavedb:00000041-a-1"], indirect=True
)
def test_for_mane_select(obtain_transcripts):
    """
    Tests for transcript selections.
    Compares computed values of RefSeq Protein and Nucleotide ID with those obtained by the notebooks
    Also checks for status of Transcripts found
    """
    # Checking for RefSeq protien ID
    computed_prot = obtain_transcripts["RefSeq_prot"]
    assert computed_prot == "NP_938033.1"

    # Checking for RefSeq Nucleotide ID
    computed_nuc = obtain_transcripts["RefSeq_nuc"]
    assert computed_nuc == "NM_198291.3"

    assert obtain_transcripts["status"] == "MANE Select"


@pytest.mark.parametrize(
    "obtain_transcripts", ["urn:mavedb:00000091-a-1"], indirect=True
)
def test_for_mane_plus_clinical(obtain_transcripts):
    """
    Tests for transcript selections.
    Compares computed values of RefSeq Protein and Nucleotide ID with those obtained by the notebooks
    Also checks for status of Transcripts found
    """
    # Checking for RefSeq protien ID
    computed_prot = obtain_transcripts["RefSeq_prot"]
    assert computed_prot == "NP_004324.2"

    # Checking for RefSeq Nucleotide ID
    computed_nuc = obtain_transcripts["RefSeq_nuc"]
    assert computed_nuc == "NM_004333.6"

    assert obtain_transcripts["status"] == "MANE Plus Clincial"


@pytest.mark.parametrize(
    "obtain_transcripts", ["urn:mavedb:00000047-c-1"], indirect=True
)
def test_for_longest_compatible(obtain_transcripts):
    """
    Tests for transcript selections.
    Compares computed values of RefSeq Protein and Nucleotide ID with those obtained by the notebooks
    Also checks for status of Transcripts found
    """
    # Checking for RefSeq protien ID
    computed_prot = obtain_transcripts["RefSeq_prot"]
    assert computed_prot == "NP_000570.1"

    # Checking for RefSeq Nucleotide ID
    computed_nuc = obtain_transcripts["RefSeq_nuc"]
    assert computed_nuc == "NM_000579.3"

    assert obtain_transcripts["status"] == "Longest Compatible"
