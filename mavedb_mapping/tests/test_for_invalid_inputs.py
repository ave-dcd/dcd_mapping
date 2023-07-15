from main import main_map
import pytest


class TestForNonHumanScoresets:
    def test_for_scoreset_if_blat_successful(self):
        scoreset = open("tests/data/urn:mavedb:00000010-a-1", "r")
        with pytest.raises(ValueError):
            main_map(scoreset)

    def test_for_scoreset_if_blat_not_successful(self):
        scoreset = open("tests/data/urn:mavedb:00000004-a-1", "r")
        with pytest.raises(ValueError):
            main_map(scoreset)


class TestForUnsuccessfulMappingScoresets:
    def test_for_no_transcripts_found(self):
        scoreset = open("tests/data/urn:mavedb:00000097-c-1", "r")
        with pytest.raises(Exception):
            main_map(scoreset)
