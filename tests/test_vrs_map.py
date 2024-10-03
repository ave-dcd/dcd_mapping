"""Test ``vrs_map.py``"""

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from cool_seq_tool.schemas import AnnotationLayer

from dcd_mapping.mavedb_data import _load_scoreset_records
from dcd_mapping.schemas import (
    AlignmentResult,
    MappedScore,
    ScoresetMetadata,
    TxSelectResult,
)
from dcd_mapping.vrs_map import vrs_map


def _assert_correct_vrs_map(
    mapping: MappedScore,
    expected_mappings_data: dict[tuple[str, AnnotationLayer], dict],
):
    key = (mapping.accession_id, mapping.annotation_layer)
    assert (
        key in expected_mappings_data
    ), "Score row/layer combination is not in expected mappings"
    expected = expected_mappings_data[key]
    assert (
        mapping.pre_mapped.id == expected["pre_mapped"]
    ), f"Mismatch in premapped ID for {key}"
    assert (
        mapping.post_mapped.id == expected["post_mapped"]
    ), f"Mismatch in postmapped ID for {key}"


@pytest.fixture()
def get_fixtures_protein(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
    transcript_results_fixture: dict[str, TxSelectResult],
):
    def _get_fixtures(urn: str):
        return (
            _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv"),
            scoreset_metadata_fixture[urn],
            align_result_fixture[urn],
            transcript_results_fixture[urn],
        )

    return _get_fixtures


@pytest.fixture()
def get_fixtures_genomic(
    fixture_data_dir: Path,
    scoreset_metadata_fixture: dict[str, ScoresetMetadata],
    align_result_fixture: dict[str, AlignmentResult],
):
    def _get_fixtures(urn: str):
        return (
            _load_scoreset_records(fixture_data_dir / f"{urn}_scores.csv"),
            scoreset_metadata_fixture[urn],
            align_result_fixture[urn],
        )

    return _get_fixtures


def test_2_a_2(
    get_fixtures_protein,
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000002-a-2"
    records, metadata, align_result, tx_result = get_fixtures_protein(urn)
    expected_mappings_data = {
        ("urn:mavedb:00000002-a-2#1", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:CPB.cbMfkw_9d-HzlhDqccrjoAs-CAvrcv2x",
            "post_mapped": "ga4gh:CPB.MHaFXDqLIiloPVIVz0MOkaQoBuFhZUCU",
        },
        ("urn:mavedb:00000002-a-2#2679", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.zBQxkMEcRjA7NXhuRTWj34sHQpsHg9y9",
            "post_mapped": "ga4gh:VA.0SXndf3nJgaSnhoJHW8ePuAmHU8vsHR2",
        },
        ("urn:mavedb:00000002-a-2#3096", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:CPB.xwtyDegestIYC_1yEPSH0QsXcsHYMDCH",
            "post_mapped": "ga4gh:CPB.JvOBJmdni3-IVkd102nv3WYdA4u-7wD1",
        },
        ("urn:mavedb:00000002-a-2#26248", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:CPB.N9DCpjE8NT04BdI8gRvLtL2UhOek0yDd",
            "post_mapped": "ga4gh:CPB.t8kUVKGpY62YLP84L5BoZLLoG-c_CiZ3",
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None
    assert len(mappings) == 4

    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        [
            "DVPLPAGWEMAKTSSGQRYFLNHIDQTTTWQDPR",
            [{"namespace": "ga4gh", "alias": "SQ.-1zvs4OMc7npphYxRnz-0lO69ueqop8R"}],
        ],
        [
            "GACGTTCCACTGCCGGCTGGTTGGGAAATGGCTAAAACTAGTTCTGGTCAGCGTTACTTCCTGAACCACATCGACCAGACCACCACGTGGCAGGACCCGCGT",
            [{"namespace": "ga4gh", "alias": "SQ.CR4vhHED4FbR2mgu0k14tfJ0ldnMoi2A"}],
        ],
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_41_a_1(
    get_fixtures_protein,
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000041-a-1"
    records, metadata, align_result, tx_result = get_fixtures_protein(urn)

    expected_mappings_data = {
        ("urn:mavedb:00000041-a-1#548", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.r6eDcCrZ6ENktOQOadid32LgOOri57WU",
            "post_mapped": "ga4gh:VA.e6KSHAit_F6PNSIGLrkhX15SbO37nBwo",
        },
        ("urn:mavedb:00000041-a-1#50", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.0IGoEHQa1DtpRiJ_M-W0JMkk1Poqbku7",
            "post_mapped": "ga4gh:VA.K29Syz7w2D-w6daR_Km0322GPOBJ6NMD",
        },
        ("urn:mavedb:00000041-a-1#51", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.9COhg0lOjgBO7IvAOSxu_nsG7ZGoZWe7",
            "post_mapped": "ga4gh:VA.c-8m280hhcSxHwz2X44g_6GOJ_LHgfWB",
        },
        ("urn:mavedb:00000041-a-1#977", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.ADQgLQLsUzwNsm-NAsev1RHAQfTS-OsB",
            "post_mapped": "ga4gh:VA.OW_YgI3iubDZhImByYN4Ukb8st_FhNUc",
        },
        ("urn:mavedb:00000041-a-1#52", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.BTiZlgO5clndCeQLggE47uMsO5cJOW2W",
            "post_mapped": "ga4gh:VA.Q3uQCBSe9xvAf0sZ9RW1VCZ4aG1wfNH1",
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
    get_fixtures_genomic,
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000099-a-1"
    records, metadata, align_result = get_fixtures_genomic(urn)

    expected_mappings_data = {
        ("urn:mavedb:00000099-a-1#8", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.PlDkRFWzI0iqCkN-LcS-fWI22kn3YZUv",
            "post_mapped": "ga4gh:VA.mcz8924gSgojiH6DIZleiDUhWZZLNGj9",
        },
        ("urn:mavedb:00000099-a-1#96", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.zx0sJRB31AojS7O_OESTwrTmVtoVrfuE",
            "post_mapped": "ga4gh:VA.ONGOnA8T9qq5kC4mTGXOOeO8JZWOJJjv",
        },
        ("urn:mavedb:00000099-a-1#194", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.-xZTr5Y325jiEWR0oMuxP-ALmjoXJFuR",
            "post_mapped": "ga4gh:VA.iT_1RID-2fnNfOo910x3cvRR454xsXYf",
        },
        ("urn:mavedb:00000099-a-1#211", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.CZL9iNU04WzGyypEHyq2kkk-7NDGgAZb",
            "post_mapped": "ga4gh:VA.lPhWnHAhGQ8zNeVXMAoM-XvB6wKLacEc",
        },
    }
    mappings = vrs_map(metadata, align_result, records, transcript=None)
    assert mappings is not None
    assert len(mappings) == 4

    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "ATGAATGGCACAGAAGGCCCTAACTTCTACGTGCCCTTCTCCAATGCGACGGGTGTGGTACGCAGCCCCTTCGAGTACCCACAGTACTACCTGGCTGAGCCATGGCAGTTCTCCATGCTGGCCGCCTACATGTTTCTGCTGATCGTGCTGGGCTTCCCCATCAACTTCCTCACGCTCTACGTCACCGTCCAGCACAAGAAGCTGCGCACGCCTCTCAACTACATCCTGCTCAACCTAGCCGTGGCTGACCTCTTCATGGTCCTAGGTGGCTTCACCAGCACCCTCTACACCTCTCTGCATGGATACTTCGTCTTCGGGCCCACAGGATGCAATTTGGAGGGCTTCTTTGCCACCCTGGGCGGTGAAATTGCCCTGTGGTCCTTGGTGGTCCTGGCCATCGAGCGGTACGTGGTGGTGTGTAAGCCCATGAGCAACTTCCGCTTCGGGGAGAACCATGCCATCATGGGCGTTGCCTTCACCTGGGTCATGGCGCTGGCCTGCGCCGCACCCCCACTCGCCGGCTGGTCCAGGTACATCCCCGAGGGCCTGCAGTGCTCGTGTGGAATCGACTACTACACGCTCAAGCCGGAGGTCAACAACGAGTCTTTTGTCATCTACATGTTCGTGGTCCACTTCACCATCCCCATGATTATCATCTTTTTCTGCTATGGGCAGCTCGTCTTCACCGTCAAGGAGGCCGCTGCCCAGCAGCAGGAGTCAGCCACCACACAGAAGGCAGAGAAGGAGGTCACCCGCATGGTCATCATCATGGTCATCGCTTTCCTGATCTGCTGGGTGCCCTACGCCAGCGTGGCATTCTACATCTTCACCCACCAGGGCTCCAACTTCGGTCCCATCTTCATGACCATCCCAGCGTTCTTTGCCAAGAGCGCCGCCATCTACAACCCTGTCATCTATATCATGATGAACAAGCAGTTCCGGAACTGCATGCTCACCACCATCTGCTGCGGCAAGAACCCACTGGGTGACGATGAGGCCTCTGCTACCGTGTCCAAGACGGAGACGAGCCAGGTGGCCCCGGCCTAA",
            [{"namespace": "ga4gh", "alias": "SQ.RtClQI6dD3uj5T-DyOr9P9_-sYR2aHJX"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_103_c_1(
    get_fixtures_protein,
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000103-c-1"
    records, metadata, align_result, tx_result = get_fixtures_protein(urn)

    expected_mappings_data = {
        ("urn:mavedb:00000103-c-1#376", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.P1WtiqLzNEm5hcxt7RFAjLBxuePK349I",
            "post_mapped": "ga4gh:VA.P1WtiqLzNEm5hcxt7RFAjLBxuePK349I",
        },
        ("urn:mavedb:00000103-c-1#55", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.TUF_eFev40uLySF3P4mrYmJEWFfjBurb",
            "post_mapped": "ga4gh:VA.TUF_eFev40uLySF3P4mrYmJEWFfjBurb",
        },
        ("urn:mavedb:00000103-c-1#2548", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.PaHvFeQptrtSZp4hXoFRWf29gkL-sppK",
            "post_mapped": "ga4gh:VA.PaHvFeQptrtSZp4hXoFRWf29gkL-sppK",
        },
        ("urn:mavedb:00000103-c-1#6810", AnnotationLayer.PROTEIN): {
            "pre_mapped": "ga4gh:VA.FhEiayb_9YGZm5tuhaWPZ65J9lfhAu8g",
            "post_mapped": "ga4gh:VA.FhEiayb_9YGZm5tuhaWPZ65J9lfhAu8g",
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=tx_result)
    assert mappings is not None
    assert len(mappings) == 4
    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIKKISPFEHQTYCQRTLREIKILLRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQHLSNDHICYFLYQILRGLKYIHSANVLHRDLKPSNLLLNTTCDLKICDFGLARVADPDHDHTGFLTEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNRPIFPGKHYLDQLNHILGILGSPSQEDLNCIINLKARNYLLSLPHKNKVPWNRLFPNADSKALDLLDKMLTFNPHKRIEVEQALAHPYLEQYYDPSDEPIAEAPFKFDMELDDLPKEKLKELIFEETARFQPGYRS",
            [{"namespace": "ga4gh", "alias": "SQ.N-m1tI22kffhKfdRZK8wCOR3QfI-1lfr"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)


def test_1_b_2(
    get_fixtures_genomic,
    mock_seqrepo_access: MagicMock,
):
    urn = "urn:mavedb:00000001-b-2"
    records, metadata, align_result = get_fixtures_genomic(urn)

    mappings = vrs_map(metadata, align_result, records)
    assert mappings is not None

    expected_mappings_data = {
        ("urn:mavedb:00000001-b-2#444", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.b8AUuR4V1jtBnur0N3uj9PApq3E9QOWX",
            "post_mapped": "ga4gh:VA.dZH6PLqHSA5sYej27eUKZ4uZ-pvg0_Q5",
        },
        ("urn:mavedb:00000001-b-2#57", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.C96gWbNLz1q7nKp3voKnd_HGDtyxdpYo",
            "post_mapped": "ga4gh:VA.EgipN3KmvY9ctrFEsUb2TdBuu48aKX4K",
        },
        ("urn:mavedb:00000001-b-2#2311", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.dh-HALPan94lHNjVY7hMWIrKvshx0sHb",
            "post_mapped": "ga4gh:VA.dY_6FNTOFoK1HvXlrqJDY-ewtdtOscnA",
        },
        ("urn:mavedb:00000001-b-2#2312", AnnotationLayer.GENOMIC): {
            "pre_mapped": "ga4gh:VA.FmoBwbceLJMI6D4U7Bt1qrBh66_oYAui",
            "post_mapped": "ga4gh:VA.H5m4fvz0i6fjkyTYkMsbQkYXmQK5dZUW",
        },
    }

    mappings = vrs_map(metadata, align_result, records, transcript=None)
    assert mappings is not None
    assert len(mappings) == 4
    for m in mappings:
        _assert_correct_vrs_map(m, expected_mappings_data)

    store_calls = [
        (
            "ATGTCTGACCAGGAGGCAAAACCTTCAACTGAGGACTTGGGGGATAAGAAGGAAGGTGAATATATTAAACTCAAAGTCATTGGACAGGATAGCAGTGAGATTCACTTCAAAGTGAAAATGACAACACATCTCAAGAAACTCAAAGAATCATACTGTCAAAGACAGGGTGTTCCAATGAATTCACTCAGGTTTCTCTTTGAGGGTCAGAGAATTGCTGATAATCATACTCCAAAAGAACTGGGAATGGAGGAAGAAGATGTGATTGAAGTTTATCAGGAACAAACGGGGGGTCATTCAACAGTTTAG",
            [{"namespace": "ga4gh", "alias": "SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk"}],
        ),
        (
            "ATGTCTGACCAGGAGGCAAAACCTTCAACTGAGGACTTGGGGGATAAGAAGGAAGGTGAATATATTAAACTCAAAGTCATTGGACAGGATAGCAGTGAGATTCACTTCAAAGTGAAAATGACAACACATCTCAAGAAACTCAAAGAATCATACTGTCAAAGACAGGGTGTTCCAATGAATTCACTCAGGTTTCTCTTTGAGGGTCAGAGAATTGCTGATAATCATACTCCAAAAGAACTGGGAATGGAGGAAGAAGATGTGATTGAAGTTTATCAGGAACAAACGGGGGGTCATTCAACAGTTTAG",
            [{"namespace": "ga4gh", "alias": "SQ.i1KiGldkfULl1XcEI-XBwhiM7x3PK5Xk"}],
        ),
    ]
    for call in store_calls:
        mock_seqrepo_access.sr.store.assert_any_call(*call)
    assert len(store_calls) == len(mock_seqrepo_access.sr.store.call_args_list)
