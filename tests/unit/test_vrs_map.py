"""Test ``vrs_map.py``

Todo:
----
* Sample a second mapping row for each test case (also, order shouldn't matter)
* Add a few more test cases

"""
from pathlib import Path
from typing import Dict

from dcd_mapping.resources import _load_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    TxSelectResult,
    VrsMapping,
)
from dcd_mapping.vrs_map import vrs_map


def test_vrs_map_93_a_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
):
    """Test VRS mapping for urn:mavedb:00000093-a-1"""
    urn = "urn:mavedb:00000093-a-1"
    scores_csv = fixture_data_dir / f"{urn}_scores.csv"
    vrs_mapping_result = vrs_map(
        scoreset_metadata_fixture[urn],
        align_result_fixture[urn],
        transcript_results_fixture[urn],
        _load_scoreset_records(scores_csv),
    )
    assert vrs_mapping_result, "Mapping result shouldn't be None"
    row: VrsMapping = vrs_mapping_result.variations[0]
    json_pre, json_post = row.output_vrs_1_3_variations()
    assert json_pre["id"] == "ga4gh:VA.XXXXXX"  # TODO
    assert json_post["id"] == "ga4gh:VA.XXXXXX"  # TODO


# def test_vrs_map_src(
#     scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
#     align_result_fixture: Dict[str, AlignmentResult],
#     transcript_results_fixture: Dict[str, TxSelectResult],
# ):
#     """Test VRS mapping for urn:mavedb:00000041-a-1"""
#     urn = "urn:mavedb:00000041-a-1"
#     vrs_mapping_result = vrs_map(
#         scoreset_metadata_fixture[urn],
#         align_result_fixture[urn],
#         transcript_results_fixture[urn],
#         get_scoreset_records(urn),
#     )
#     assert vrs_mapping_result is not None
#     mapping = vrs_mapping_result.variations[0]
#     assert mapping.mavedb_id == "urn:mavedb:00000041-a-1#548"
#     assert mapping.score == pytest.approx(1.513315421)
#
#     assert isinstance(mapping.pre_mapped, Allele)
#     assert mapping.pre_mapped.state.sequence == "T"
#     assert (
#         mapping.pre_mapped.location.sequenceReference.refgetAccession
#         == "ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB"
#     )
#     assert mapping.pre_mapped.location.start == 14
#     assert mapping.pre_mapped.location.end == 15
#
#     assert isinstance(mapping.post_mapped, Allele)
#     assert mapping.post_mapped.state.sequence == "T"
#     assert (
#         mapping.post_mapped.location.sequenceReference.refgetAccession
#         == "ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye"
#     )
#     assert mapping.post_mapped.location.start == 283
#     assert mapping.post_mapped.location.end == 284
#     assert mapping.post_mapped.expressions[0]["syntax"] == "hgvs.p"
#     assert mapping.post_mapped.expressions[0]["value"] == "NNP_005408.1:p.Val284Thr"
#     assert mapping.post_mapped.expressions[0]["syntax_version"] is None
#
#     assert mapping.mapped_transcript is not None
#     assert mapping.mapped_transcript.refseq_nuc == "NM_198291.3"
#     assert (
#         mapping.mapped_transcript.transcript_priority == TranscriptPriority.MANE_SELECT
#     )


# def test_vrs_map_hbb(
#     scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
#     align_result_fixture: Dict[str, AlignmentResult],
#     transcript_results_fixture: Dict[str, TxSelectResult],
# ):
#     """Test VRS mapping for urn:mavedb:00000018-a-1"""
#     urn = "urn:mavedb_id:00000018-a-1"
#     vrs_mapping_result = vrs_map(
#         scoreset_metadata_fixture[urn],
#         align_result_fixture[urn],
#         transcript_results_fixture[urn],
#         get_scoreset_records(urn),
#     )
#     assert vrs_mapping_result is not None
#     mapping = vrs_mapping_result.variations[0]
#     assert mapping.mavedb_id == "urn:mavedb:00000018-a-1#583"
#     assert mapping.score == pytest.approx(-0.17)
#
#     assert isinstance(mapping.pre_mapped, Allele)
#     assert mapping.pre_mapped.state.sequence == "G"
#     assert (
#         mapping.pre_mapped.location.sequenceReference.refgetAccession
#         == "ga4gh:SQ.jUOcLPDjSqWFEo9kSOG8ITe1dr9QK3h6"
#     )
#     assert mapping.pre_mapped.location.start == 186
#     assert mapping.pre_mapped.location.end == 187
#
#     assert isinstance(mapping.post_mapped, Allele)
#     assert mapping.post_mapped.state.sequence == "G"
#     assert (
#         mapping.post_mapped.location.sequenceReference.refgetAccession
#         == "ga4gh:SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1"
#     )
#     assert mapping.post_mapped.location.start == 5227207
#     assert mapping.post_mapped.location.end == 5227208
#     assert mapping.post_mapped.expressions[0]["syntax"] == "hgvs.g"
#     assert mapping.post_mapped.expressions[0]["value"] == "NC_000011.10:g.5227208T>G"
#     assert mapping.post_mapped.expressions[0]["syntax_version"] is None
#
#     assert mapping.mapped_transcript is not None
#     assert mapping.mapped_transcript.refseq_nuc == "NM_198291.3"
#     assert (
#         mapping.mapped_transcript.transcript_priority == TranscriptPriority.MANE_SELECT
#     )
