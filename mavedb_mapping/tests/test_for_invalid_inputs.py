from main import main_map
import pytest

"""Tests for scoresets that raise errors/unsuccessful mappings."""


class TestForNonHumanScoresets:

    """Test Class to test for non human scoresets that are expected to raise errors"""

    def test_for_scoreset_if_blat_successful(self):
        """Test to check for non human scoreset where BLAT Alignment is successful, and gives low coverage.
        The mapping raises an Exception"""

        scoreset = open("tests/data/urn:mavedb:00000010-a-1", "r")
        with pytest.raises(Exception):
            main_map(scoreset)

    def test_for_scoreset_if_blat_not_successful(self):
        """Test to check for non human scoreset where BLAT Alignment is not successful, and gives no output.
        The mapping raises a ValueError"""

        scoreset = open("tests/data/urn:mavedb:00000004-a-1", "r")
        with pytest.raises(ValueError):
            main_map(scoreset)


class TestForUnsuccessfulMappingScoresets:

    """Test class for human scoresets with unsuccessful mapping"""

    def test_for_no_transcripts_found(self):
        """Testing for a scoreset for which transcripts could not be found"""
        # TODO: create test for 53
        scoreset = open("tests/data/urn:mavedb:00000097-c-1", "r")
        with pytest.raises(Exception):
            main_map(scoreset)
