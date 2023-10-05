"""Test ``align`` module."""
import json
from pathlib import Path

import pytest
from mavemap.align import align
from mavemap.schemas import ScoresetMetadata


@pytest.fixture(scope="module")
def scoreset_metadata_fixture():
    """Provide scoreset metadata fixtures."""
    fixture_file = (
        Path(__file__).parents[1].resolve() / "fixtures" / "scoreset_metadata.json"
    )
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for d in data["scoreset_metadata"]:
        formatted_data = ScoresetMetadata(**d)
        results[formatted_data.target_gene_name] = formatted_data
    return results


def test_align(scoreset_metadata_fixture):
    """Test ``align()`` method.

    We should be able to run this without an available BLAT binary, so we'll mock the
    ``_run_blat_command()`` call.
    """
    # urn:mavedb:00000041-a-1
    # Src catalytic domain
    scoreset_metadata = scoreset_metadata_fixture["Src catalytic domain"]
    align_result = align(scoreset_metadata)
    assert align_result
    assert align_result.chrom == "20"
    assert align_result.strand == 1
    assert align_result.coverage == pytest.approx(100.0)
    assert align_result.query_range.start == 0
    assert align_result.query_range.end == 750
    query_subranges = [
        [0, 52],
        [52, 232],
        [232, 309],
        [309, 463],
        [463, 595],
        [595, 750],
    ]
    for actual, expected in zip(align_result.query_subranges, query_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
    assert align_result.hit_range.start == 37397802
    assert align_result.hit_range.end == 37403325
    hit_subranges = [
        [37397802, 37397854],
        [37400114, 37400294],
        [37401601, 37401678],
        [37402434, 37402588],
        [37402748, 37402880],
        [37403170, 37403325],
    ]
    for actual, expected in zip(align_result.hit_subranges, hit_subranges):
        assert actual.start == expected[0]
        assert actual.end == expected[1]
