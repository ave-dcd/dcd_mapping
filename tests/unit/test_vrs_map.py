"""Test ``vrs_map.py``"""
from pathlib import Path
from typing import Dict, Optional
from unittest.mock import MagicMock

import pytest

from dcd_mapping.mavedb_data import _load_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    TxSelectResult,
    VrsObject1_x,
)
from dcd_mapping.vrs_map import vrs_map


def _assert_correct_vrs_map(
    mapping: VrsObject1_x, expected_mappings_data: Dict[str, Dict]
):
    assert (
        mapping.mavedb_id in expected_mappings_data
    ), "Score row is in expected mappings"
    expected = expected_mappings_data[mapping.mavedb_id]
    assert mapping.pre_mapped_variants["id"] == expected["pre_mapped"]
    assert mapping.post_mapped_variants["id"] == expected["post_mapped"]


@pytest.fixture()
def mock_seqrepo(mocker: MagicMock):
    """Mock SeqRepo instance.

    So far, it seems like the DataProxy `get_sequence` method is all that needs mocking;
    to add new test cases, throw some print/breakpoint statements into the original
    methods and capture results there.
    """
    mock_seqrepo_instance = mocker.MagicMock()

    def _get_sequence(
        identifier: str, start: Optional[str] = None, end: Optional[str] = None
    ) -> str:
        calls = {
            # 41-a-1
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 14, 15): "V",
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 103, 104): "I",
            ("ga4gh:SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB", 149, 150): "Y",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 14, 15): "R",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 103, 104): "S",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 149, 150): "E",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 283, 284): "V",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 372, 373): "I",
            ("ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye", 418, 419): "Y",
            # 103-c-1
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 336, 337): "D",
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 352, 353): "R",
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 220, 221): "M",
            ("ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr", 1, 2): "A",
            # 1-b-2 (WIP)
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 75, 76): "T",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 225, 228): "ACT",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202210743, 202210746): "AGT",
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 7, 8): "P",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 21, 24): "CCT",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 21, 23): "CC",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202220094, 202220097): "AGG",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202220095, 202220097): "GG",
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 68, 69): "Q",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 204, 207): "CAG",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 205, 207): "AG",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202210764, 202210767): "CTG",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202210764, 202210766): "CT",
            ("ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip", 81, 82): "M",
            ("ga4gh:SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk", 245, 246): "G",
            ("ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g", 202207312, 202207313): "C",
        }
        return calls[(identifier, start, end)]

    def _translate_sequence_identifier(
        identifier: str, namespace: Optional[str] = None
    ) -> str:
        calls = {
            ("refseq:NP_938033.1", "ga4gh"): [
                "ga4gh:SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye"
            ],
            ("refseq:NP_003343.1", "ga4gh"): [
                "ga4gh:SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip"
            ],
            ("refseq:NC_000002.12", "ga4gh"): [
                "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g"
            ],
            ("refseq:NP_002736.3", "ga4gh"): [
                "ga4gh:SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr"
            ],
        }
        return calls[(identifier, namespace)]

    def _derive_refget_accession(ac: str) -> str:
        calls = {
            "NP_938033.1": "SQ.uJDQo_HaTNFL2-0-6K5dVzVcweigexye",
            "NP_002736.3": "SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr",
        }
        return calls[ac]

    mock_seqrepo_access = mocker.MagicMock()
    mock_seqrepo_access.get_sequence.side_effect = _get_sequence
    mock_seqrepo_access.translate_sequence_identifier.side_effect = (
        _translate_sequence_identifier
    )
    mock_seqrepo_access.derive_refget_accession.side_effect = _derive_refget_accession
    mock_seqrepo_access.sr = mock_seqrepo_instance
    mocker.patch("dcd_mapping.vrs_map.get_seqrepo", return_value=mock_seqrepo_access)
    mocker.patch("dcd_mapping.lookup.get_seqrepo", return_value=mock_seqrepo_access)
    return mock_seqrepo_access


def test_41_a_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
    mock_seqrepo: MagicMock,
):
    urn = "urn:mavedb:00000041-a-1"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_result = transcript_results_fixture[urn]

    expected_mappings_data = {
        "urn:mavedb:00000041-a-1#548": {
            "pre_mapped": "ga4gh:VA.NJgaCF0JPFERdw9Y7fW4bXkaP1tSa5fv",
            "post_mapped": "ga4gh:VA.csCB31gWoiiD38TlR35dZnyAI156YWgW",
        },
        "urn:mavedb:00000041-a-1#50": {
            "pre_mapped": "ga4gh:VA.RfNyaPcZg8o9f3YK0qa4Vu3LwBEQdmVf",
            "post_mapped": "ga4gh:VA.QP9KLxvN6_b7sWeY7L8FBs5XKMxsCiLE",
        },
        "urn:mavedb:00000041-a-1#51": {
            "pre_mapped": "ga4gh:VA.sZNa3SNPlv_gU2JSiH7Q03nNT7oFy1NX",
            "post_mapped": "ga4gh:VA.2kmNV4T4Bp_UAnI002QOfrzd_yqb21vs",
        },
        "urn:mavedb:00000041-a-1#977": {
            "pre_mapped": "ga4gh:VA.h7QtWm0WzlOr0zbB9y4tJOIEUet6VLcB",
            "post_mapped": "ga4gh:VA.bakDMkeUFIa46_HAwXb7gUjhywggVNIN",
        },
        "urn:mavedb:00000041-a-1#52": {
            "pre_mapped": "ga4gh:VA.gl5xiWNmwUfEMZe5Aub15HIiztUaizay",
            "post_mapped": "ga4gh:VA.3Pp5-tRnYkmm8f6qxk06GvTpn81DqiQV",
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None

    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        [
            "LRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGETGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKEPEERPTFEYLQAFL",
            [{"namespace": "ga4gh", "alias": "SQ.PyX9IDu95_tYLg1Jz9JpW5xpQkwn6bpB"}],
        ],
        [
            "CTGCGGCTGGAGGTCAAGCTGGGCCAGGGCTGCTTTGGCGAGGTGTGGATGGGGACCTGGAACGGTACCACCAGGGTGGCCATCAAAACCCTGAAGCCTGGCACGATGTCTCCAGAGGCCTTCCTGCAGGAGGCCCAGGTCATGAAGAAGCTGAGGCATGAGAAGCTGGTGCAGTTGTATGCTGTGGTTTCAGAGGAGCCCATTTACATCGTCACGGAGTACATGAGCAAGGGGAGTTTGCTGGACTTTCTCAAGGGGGAGACAGGCAAGTACCTGCGGCTGCCTCAGCTGGTGGACATGGCTGCTCAGATCGCCTCAGGCATGGCGTACGTGGAGCGGATGAACTACGTCCACCGGGACCTTCGTGCAGCCAACATCCTGGTGGGAGAGAACCTGGTGTGCAAAGTGGCCGACTTTGGGCTGGCTCGGCTCATTGAAGACAATGAGTACACGGCGCGGCAAGGTGCCAAATTCCCCATCAAGTGGACGGCTCCAGAAGCTGCCCTCTATGGCCGCTTCACCATCAAGTCGGACGTGTGGTCCTTCGGGATCCTGCTGACTGAGCTCACCACAAAGGGACGGGTGCCCTACCCTGGGATGGTGAACCGCGAGGTGCTGGACCAGGTGGAGCGGGGCTACCGGATGCCCTGCCCGCCGGAGTGTCCCGAGTCCCTGCACGACCTCATGTGCCAGTGCTGGCGGAAGGAGCCTGAGGAGCGGCCCACCTTCGAGTACCTGCAGGCCTTCCTG",
            [{"namespace": "ga4gh", "alias": "SQ.SBlyhKAQ5u14huYar7I7UmXZk4eRPexx"}],
        ],
    ]
    for call in store_calls:
        mock_seqrepo.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo.sr.store.call_args_list)


def test_103_c_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
    transcript_results_fixture: Dict[str, TxSelectResult],
    mock_seqrepo: MagicMock,
):
    urn = "urn:mavedb:00000103-c-1"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_result = transcript_results_fixture[urn]

    expected_mappings_data = {
        "urn:mavedb:00000103-c-1#376": {
            "pre_mapped": "ga4gh:VA.KPfJzFOpDl49yIPqqeF07BsnCsOPbf5Z",
            "post_mapped": "ga4gh:VA.KPfJzFOpDl49yIPqqeF07BsnCsOPbf5Z",
        },
        "urn:mavedb:00000103-c-1#55": {
            "pre_mapped": "ga4gh:VA.AVg1O4zA7DP71z-QZD7-5864A52Cp4RO",
            "post_mapped": "ga4gh:VA.AVg1O4zA7DP71z-QZD7-5864A52Cp4RO",
        },
        "urn:mavedb:00000103-c-1#2548": {
            "pre_mapped": "ga4gh:VA.4dd3ml7ZuqyoMhwJhMOWi2n1729MtE-6",
            "post_mapped": "ga4gh:VA.4dd3ml7ZuqyoMhwJhMOWi2n1729MtE-6",
        },
        "urn:mavedb:00000103-c-1#6810": {
            "pre_mapped": "ga4gh:VA.H_7xCjw--slTAAKFwP42WCUo1tITw39d",
            "post_mapped": "ga4gh:VA.H_7xCjw--slTAAKFwP42WCUo1tITw39d",
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None

    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIKKISPFEHQTYCQRTLREIKILLRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQHLSNDHICYFLYQILRGLKYIHSANVLHRDLKPSNLLLNTTCDLKICDFGLARVADPDHDHTGFLTEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNRPIFPGKHYLDQLNHILGILGSPSQEDLNCIINLKARNYLLSLPHKNKVPWNRLFPNADSKALDLLDKMLTFNPHKRIEVEQALAHPYLEQYYDPSDEPIAEAPFKFDMELDDLPKEKLKELIFEETARFQPGYRS",
            [{"namespace": "ga4gh", "alias": "SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr"}],
        )
    ]
    for call in store_calls:
        mock_seqrepo.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo.sr.store.call_args_list)


# def test_1_b_2(
#     fixture_data_dir: Path,
#     scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
#     align_result_fixture: Dict[str, AlignmentResult],
#     transcript_results_fixture: Dict[str, TxSelectResult],
#     mock_seqrepo: MagicMock,
# ):
#     urn = "urn:mavedb:00000001-b-2"
#     records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
#     metadata = scoreset_metadata_fixture[urn]
#     align_result = align_result_fixture[urn]
#     tx_result = transcript_results_fixture[urn]
#
#     mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
#     assert mappings is not None
#
#     expected_mappings_data = {
#         "urn:mavedb:00000001-b-2#444": {
#             "pre_mapped": "ga4gh:VA.1XfnARSJMWLoHTQDJU-m0ZYLUsj-qRNM",
#             "post_mapped": "ga4gh:VA.ifSwfAlXaZWIcTQaXrZv7LCsa6sywRvw",
#         },
#         "urn:mavedb:00000001-b-2#57": {
#             "pre_mapped": "ga4gh:VA.-zYenjM1Wsuu5Ia06nn856cn4JKEjMwR",
#             "post_mapped": "ga4gh:VA.LODyOWFdnBsGn7dciC-MCSLHCNtoyjaf",
#         },
#         "urn:mavedb:00000001-b-2#2311": {
#             "pre_mapped": "ga4gh:VA.kdSdynoQvgoau0GqsUzhduF7riQVvx5-",
#             "post_mapped": "ga4gh:VA.LdKB-BHsNMueerA0u4ngGU1oxK2oBDHs",
#         },
#         "urn:mavedb:00000001-b-2#2312": {
#             "pre_mapped": "ga4gh:VA.pkEhKe6fG7eWaUTA7DJG2v-zggeP5N1m",
#             "post_mapped": "ga4gh:VA.R4iEL0X_2Mr4o4cmuLE4AW_UrTotjZ8M",
#         },
#     }
#
#     mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
#     assert mappings is not None
#
#     def _assert_correct_vrs_map(mapping: VrsObject1_x):
#         assert (
#             mapping.mavedb_id in expected_mappings_data
#         ), "Score row is in expected mappings"
#         expected = expected_mappings_data[mapping.mavedb_id]
#         assert mapping.pre_mapped_variants["id"] == expected["pre_mapped"]
#         assert mapping.post_mapped_variants["id"] == expected["post_mapped"]
#
#     for m in mappings[:2]:
#         _assert_correct_vrs_map(m)
#
#     store_calls = [
#         (
#             "ATGTCTGACCAGGAGGCAAAACCTTCAACTGAGGACTTGGGGGATAAGAAGGAAGGTGAATATATTAAACTCAAAGTCATTGGACAGGATAGCAGTGAGATTCACTTCAAAGTGAAAATGACAACACATCTCAAGAAACTCAAAGAATCATACTGTCAAAGACAGGGTGTTCCAATGAATTCACTCAGGTTTCTCTTTGAGGGTCAGAGAATTGCTGATAATCATACTCCAAAAGAACTGGGAATGGAGGAAGAAGATGTGATTGAAGTTTATCAGGAACAAACGGGGGGTCATTCAACAGTTTAG",
#             [{"namespace": "ga4gh", "alias": "SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk"}],
#         ),
#         (
#             "sdlkfsdlksdflk",  # TODO -- update once this test is passing above
#             [{"namespace": "ga4gh", "alias": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g"}],
#         ),
#     ]
#     for call in store_calls:
#         mock_seqrepo.sr.store.assert_any_call(*call)
