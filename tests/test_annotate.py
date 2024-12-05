"""Test annotation module."""

from unittest.mock import MagicMock

import pytest

from dcd_mapping.annotate import annotate
from dcd_mapping.schemas import (
    AlignmentResult,
    MappedScore,
    ScoresetMetadata,
    TxSelectResult,
)


@pytest.fixture()
def get_fixtures(
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    transcript_results_fixture: dict[str, TxSelectResult],
    mapped_scores_fixture: dict[str, list[MappedScore]],
    align_result_fixture: dict[str, AlignmentResult],
):
    def _get_fixtures(urn: str):
        return (
            mapped_scores_fixture[urn],
            transcript_results_fixture[urn],
            scoreset_metadata_fixture[urn],
            align_result_fixture[urn],
        )

    return _get_fixtures


def test_2_a_2(get_fixtures, mock_seqrepo_access: MagicMock):  # noqa: ARG001
    urn = "urn:mavedb:00000002-a-2"
    mapped_scores, tx_results, metadata, align_result = get_fixtures(urn)

    annotate_result = annotate(mapped_scores, tx_results, metadata, align_result)

    expected_list = [
        {
            "pre_1_3_id": "ga4gh:VH.jb2M9uN8u-YZ9x9mTbuinFlpQjXqDLF2",
            "post_1_3_id": "ga4gh:VH.ufHLQtHCsVphPAt1aD5XKO43JS8BQkgD",
        },
        {
            "pre_1_3_id": "ga4gh:VA.lPVV_PjqY_htCoCjblDHH9kX2S6V9W-R",
            "post_1_3_id": "ga4gh:VA.e73a5KuabQy1-2BatNKde1mTm0G3Yc2l",
        },
        {
            "pre_1_3_id": "ga4gh:VH.RE7uz4k7WxhGfH6r41388YcXCeYE7Qcc",
            "post_1_3_id": "ga4gh:VH.nRXhn2igvfh8a0PPg0QtgGGFdWQSFIky",
        },
        {
            "pre_1_3_id": "ga4gh:VH.6RgmRfEehvXZchrSKSOx8oeUatodglNZ",
            "post_1_3_id": "ga4gh:VH.fKALOa7vJHwMJUgwVV1SL9Uiq02B0be9",
        },
    ]

    for actual, expected in zip(annotate_result, expected_list, strict=False):
        pre_1_3_id = next(
            e for e in actual.pre_mapped.extensions if e.name == "vrs_v1.3_id"
        ).value
        assert pre_1_3_id == expected["pre_1_3_id"]
        post_1_3_id = next(
            e for e in actual.post_mapped.extensions if e.name == "vrs_v1.3_id"
        ).value
        assert post_1_3_id == expected["post_1_3_id"]
