"""Test ``vrs_map.py``"""
from pathlib import Path
from unittest.mock import MagicMock

from cool_seq_tool.schemas import AnnotationLayer

from dcd_mapping.mavedb_data import _load_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    TxSelectResult,
    VrsObject1_x,
)
from dcd_mapping.vrs_map import vrs_map


def _assert_correct_vrs_map(
    mapping: VrsObject1_x, expected_mappings_data: dict[str, dict]
):
    assert (
        mapping.mavedb_id in expected_mappings_data
    ), "Score row is in expected mappings"
    assert (
        mapping.layer in expected_mappings_data[mapping.mavedb_id]
    ), "Includes expected mapping layer for score row"
    expected = expected_mappings_data[mapping.mavedb_id][mapping.layer]
    assert mapping.pre_mapped_variants["id"] == expected["pre_mapped"]
    assert mapping.post_mapped_variants["id"] == expected["post_mapped"]


def test_41_a_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000041-a-1"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_result = transcript_results_fixture[urn]

    expected_mappings_data = {
        "urn:mavedb:00000041-a-1#548": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.NJgaCF0JPFERdw9Y7fW4bXkaP1tSa5fv",
                "post_mapped": "ga4gh:VA.csCB31gWoiiD38TlR35dZnyAI156YWgW",
            }
        },
        "urn:mavedb:00000041-a-1#50": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.RfNyaPcZg8o9f3YK0qa4Vu3LwBEQdmVf",
                "post_mapped": "ga4gh:VA.QP9KLxvN6_b7sWeY7L8FBs5XKMxsCiLE",
            }
        },
        "urn:mavedb:00000041-a-1#51": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.sZNa3SNPlv_gU2JSiH7Q03nNT7oFy1NX",
                "post_mapped": "ga4gh:VA.2kmNV4T4Bp_UAnI002QOfrzd_yqb21vs",
            }
        },
        "urn:mavedb:00000041-a-1#977": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.h7QtWm0WzlOr0zbB9y4tJOIEUet6VLcB",
                "post_mapped": "ga4gh:VA.bakDMkeUFIa46_HAwXb7gUjhywggVNIN",
            }
        },
        "urn:mavedb:00000041-a-1#52": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.gl5xiWNmwUfEMZe5Aub15HIiztUaizay",
                "post_mapped": "ga4gh:VA.3Pp5-tRnYkmm8f6qxk06GvTpn81DqiQV",
            }
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None
    assert len(mappings) == 5

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
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_99_a_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000099-a-1"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_result = transcript_results_fixture[urn]

    expected_mappings_data = {
        "urn:mavedb:00000099-a-1#8": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.ES3p2PRZ6Lhial9xnlVxV8uXUbl03ESx",
                "post_mapped": "ga4gh:VA.ES3p2PRZ6Lhial9xnlVxV8uXUbl03ESx",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.dw7KwkqdRsRL5UU55KznGYFqGCK6sWgt",
                "post_mapped": "ga4gh:VA.T8SkVZHuAF9J23etkAXG2Sz2w2yJMSO4",
            },
        },
        "urn:mavedb:00000099-a-1#96": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.bu9RmB2kvJiTFKfThhtpYam50_GrIkbx",
                "post_mapped": "ga4gh:VA.bu9RmB2kvJiTFKfThhtpYam50_GrIkbx",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.4cQNoG_u6yGDk_GMmVvbYX-8EniN6RKd",
                "post_mapped": "ga4gh:VA.4-NsddlZX70Jzi7OxSL6KYRe8T-h4rwT",
            },
        },
        "urn:mavedb:00000099-a-1#194": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.eT6tN7HdCmDQnTHuetRHtRZRHPxpO53x",
                "post_mapped": "ga4gh:VA.eT6tN7HdCmDQnTHuetRHtRZRHPxpO53x",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.gRO5ANz9hJ6y6NuyzHOzl94b39SAbbIE",
                "post_mapped": "ga4gh:VA.JSmTABIybCChEE6mmjnPuubOt9eRWj_4",
            },
        },
        "urn:mavedb:00000099-a-1#211": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.MdeuhDWe0mccL5lPkb5J3d_I39m-Zg1R",
                "post_mapped": "ga4gh:VA.MdeuhDWe0mccL5lPkb5J3d_I39m-Zg1R",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.K4bnyLQ5M3WCmli8lvAqHfqkZuk4KvwP",
                "post_mapped": "ga4gh:VA.HSfipwsg28LbqwITCawzumz_OWZYu_jM",
            },
        },
    }
    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None
    assert len(mappings) == 8  # includes protein and genomic for all 4 rows

    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA",
            [{"namespace": "ga4gh", "alias": "SQ.kntF3MvVMQFzU94f0QegL_ktmT5f2Dk9"}],
        ),
        (
            "ATGAATGGCACAGAAGGCCCTAACTTCTACGTGCCCTTCTCCAATGCGACGGGTGTGGTACGCAGCCCCTTCGAGTACCCACAGTACTACCTGGCTGAGCCATGGCAGTTCTCCATGCTGGCCGCCTACATGTTTCTGCTGATCGTGCTGGGCTTCCCCATCAACTTCCTCACGCTCTACGTCACCGTCCAGCACAAGAAGCTGCGCACGCCTCTCAACTACATCCTGCTCAACCTAGCCGTGGCTGACCTCTTCATGGTCCTAGGTGGCTTCACCAGCACCCTCTACACCTCTCTGCATGGATACTTCGTCTTCGGGCCCACAGGATGCAATTTGGAGGGCTTCTTTGCCACCCTGGGCGGTGAAATTGCCCTGTGGTCCTTGGTGGTCCTGGCCATCGAGCGGTACGTGGTGGTGTGTAAGCCCATGAGCAACTTCCGCTTCGGGGAGAACCATGCCATCATGGGCGTTGCCTTCACCTGGGTCATGGCGCTGGCCTGCGCCGCACCCCCACTCGCCGGCTGGTCCAGGTACATCCCCGAGGGCCTGCAGTGCTCGTGTGGAATCGACTACTACACGCTCAAGCCGGAGGTCAACAACGAGTCTTTTGTCATCTACATGTTCGTGGTCCACTTCACCATCCCCATGATTATCATCTTTTTCTGCTATGGGCAGCTCGTCTTCACCGTCAAGGAGGCCGCTGCCCAGCAGCAGGAGTCAGCCACCACACAGAAGGCAGAGAAGGAGGTCACCCGCATGGTCATCATCATGGTCATCGCTTTCCTGATCTGCTGGGTGCCCTACGCCAGCGTGGCATTCTACATCTTCACCCACCAGGGCTCCAACTTCGGTCCCATCTTCATGACCATCCCAGCGTTCTTTGCCAAGAGCGCCGCCATCTACAACCCTGTCATCTATATCATGATGAACAAGCAGTTCCGGAACTGCATGCTCACCACCATCTGCTGCGGCAAGAACCCACTGGGTGACGATGAGGCCTCTGCTACCGTGTCCAAGACGGAGACGAGCCAGGTGGCCCCGGCCTAA",
            [{"namespace": "ga4gh", "alias": "SQ.RtClQI6dD3uj5T-DyOr9P9_-sYR2aHJX"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_103_c_1(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000103-c-1"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_result = transcript_results_fixture[urn]

    expected_mappings_data = {
        "urn:mavedb:00000103-c-1#376": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.KPfJzFOpDl49yIPqqeF07BsnCsOPbf5Z",
                "post_mapped": "ga4gh:VA.KPfJzFOpDl49yIPqqeF07BsnCsOPbf5Z",
            }
        },
        "urn:mavedb:00000103-c-1#55": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.AVg1O4zA7DP71z-QZD7-5864A52Cp4RO",
                "post_mapped": "ga4gh:VA.AVg1O4zA7DP71z-QZD7-5864A52Cp4RO",
            }
        },
        "urn:mavedb:00000103-c-1#2548": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.4dd3ml7ZuqyoMhwJhMOWi2n1729MtE-6",
                "post_mapped": "ga4gh:VA.4dd3ml7ZuqyoMhwJhMOWi2n1729MtE-6",
            }
        },
        "urn:mavedb:00000103-c-1#6810": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.H_7xCjw--slTAAKFwP42WCUo1tITw39d",
                "post_mapped": "ga4gh:VA.H_7xCjw--slTAAKFwP42WCUo1tITw39d",
            }
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None
    assert len(mappings) == 4
    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    # makes the same `store()` call twice
    store_calls = [
        (
            "MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIKKISPFEHQTYCQRTLREIKILLRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQHLSNDHICYFLYQILRGLKYIHSANVLHRDLKPSNLLLNTTCDLKICDFGLARVADPDHDHTGFLTEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNRPIFPGKHYLDQLNHILGILGSPSQEDLNCIINLKARNYLLSLPHKNKVPWNRLFPNADSKALDLLDKMLTFNPHKRIEVEQALAHPYLEQYYDPSDEPIAEAPFKFDMELDDLPKEKLKELIFEETARFQPGYRS",
            [{"namespace": "ga4gh", "alias": "SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr"}],
        ),
        (
            "MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIKKISPFEHQTYCQRTLREIKILLRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQHLSNDHICYFLYQILRGLKYIHSANVLHRDLKPSNLLLNTTCDLKICDFGLARVADPDHDHTGFLTEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNRPIFPGKHYLDQLNHILGILGSPSQEDLNCIINLKARNYLLSLPHKNKVPWNRLFPNADSKALDLLDKMLTFNPHKRIEVEQALAHPYLEQYYDPSDEPIAEAPFKFDMELDDLPKEKLKELIFEETARFQPGYRS",
            [{"namespace": "ga4gh", "alias": "SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_1_b_2(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000001-b-2"
    records = _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv")
    metadata = scoreset_metadata_fixture[urn]
    align_result = align_result_fixture[urn]
    tx_result = transcript_results_fixture[urn]

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None

    expected_mappings_data = {
        "urn:mavedb:00000001-b-2#444": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.ojIs4GEPbxiMt5E6InnF6k6m9ix_z3SH",
                "post_mapped": "ga4gh:VA.ojIs4GEPbxiMt5E6InnF6k6m9ix_z3SH",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.1XfnARSJMWLoHTQDJU-m0ZYLUsj-qRNM",
                "post_mapped": "ga4gh:VA.ifSwfAlXaZWIcTQaXrZv7LCsa6sywRvw",
            },
        },
        "urn:mavedb:00000001-b-2#57": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.EwrANE7HrTJ3OidbdJcL4DAfiYiAJ0Zn",
                "post_mapped": "ga4gh:VA.EwrANE7HrTJ3OidbdJcL4DAfiYiAJ0Zn",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.-zYenjM1Wsuu5Ia06nn856cn4JKEjMwR",
                "post_mapped": "ga4gh:VA.LODyOWFdnBsGn7dciC-MCSLHCNtoyjaf",
            },
        },
        "urn:mavedb:00000001-b-2#2311": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.KQrLWh1WBJNOO7tEkJ7ujkkHCEGGTSCo",
                "post_mapped": "ga4gh:VA.KQrLWh1WBJNOO7tEkJ7ujkkHCEGGTSCo",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.kdSdynoQvgoau0GqsUzhduF7riQVvx5-",
                "post_mapped": "ga4gh:VA.LdKB-BHsNMueerA0u4ngGU1oxK2oBDHs",
            },
        },
        "urn:mavedb:00000001-b-2#2312": {
            AnnotationLayer.PROTEIN: {
                "pre_mapped": "ga4gh:VA.ZdHxEMecv2kdYzRJz5iaqiS5qUd05Wim",
                "post_mapped": "ga4gh:VA.ZdHxEMecv2kdYzRJz5iaqiS5qUd05Wim",
            },
            AnnotationLayer.GENOMIC: {
                "pre_mapped": "ga4gh:VA.pkEhKe6fG7eWaUTA7DJG2v-zggeP5N1m",
                "post_mapped": "ga4gh:VA.R4iEL0X_2Mr4o4cmuLE4AW_UrTotjZ8M",
            },
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None
    assert len(mappings) == 8
    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "MSDQEAKPSTEDLGDKKEGEYIKLKVIGQDSSEIHFKVKMTTHLKKLKESYCQRQGVPMNSLRFLFEGQRIADNHTPKELGMEEEDVIEVYQEQTGGHSTV",
            [{"namespace": "ga4gh", "alias": "SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip"}],
        ),
        (
            "ATGTCTGACCAGGAGGCAAAACCTTCAACTGAGGACTTGGGGGATAAGAAGGAAGGTGAATATATTAAACTCAAAGTCATTGGACAGGATAGCAGTGAGATTCACTTCAAAGTGAAAATGACAACACATCTCAAGAAACTCAAAGAATCATACTGTCAAAGACAGGGTGTTCCAATGAATTCACTCAGGTTTCTCTTTGAGGGTCAGAGAATTGCTGATAATCATACTCCAAAAGAACTGGGAATGGAGGAAGAAGATGTGATTGAAGTTTATCAGGAACAAACGGGGGGTCATTCAACAGTTTAG",
            [{"namespace": "ga4gh", "alias": "SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk"}],
        ),
        (
            "MSDQEAKPSTEDLGDKKEGEYIKLKVIGQDSSEIHFKVKMTTHLKKLKESYCQRQGVPMNSLRFLFEGQRIADNHTPKELGMEEEDVIEVYQEQTGGHSTV",
            [{"namespace": "ga4gh", "alias": "SQ.VkCzFNsbifqfq61Mud6oGmz0ID6CLIip"}],
        ),
        (
            "ATGTCTGACCAGGAGGCAAAACCTTCAACTGAGGACTTGGGGGATAAGAAGGAAGGTGAATATATTAAACTCAAAGTCATTGGACAGGATAGCAGTGAGATTCACTTCAAAGTGAAAATGACAACACATCTCAAGAAACTCAAAGAATCATACTGTCAAAGACAGGGTGTTCCAATGAATTCACTCAGGTTTCTCTTTGAGGGTCAGAGAATTGCTGATAATCATACTCCAAAAGAACTGGGAATGGAGGAAGAAGATGTGATTGAAGTTTATCAGGAACAAACGGGGGGTCATTCAACAGTTTAG",
            [{"namespace": "ga4gh", "alias": "SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)
