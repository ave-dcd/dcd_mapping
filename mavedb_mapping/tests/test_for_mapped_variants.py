# prot and non prot
import pytest
import json


def test_for_vrs_mapping_prot_coding(full_mapping):
    """Test function to compare results from sample mappings obtained from the notebooks"""

    mapped_variant, urn = full_mapping

    # Sample Mapping file
    mapped_file = open(f"tests/data/{urn}.json")
    mapped_example = json.load(mapped_file)

    # Sequence id from sample mapping file
    sample_sequence_id = mapped_example["mapped_scores"][0]["pre_mapped"]["location"][
        "sequence_id"
    ]

    # Sequence ID from mapped variant
    sequence_id = mapped_variant[urn][0]["pre_mapping"][0]["location"]["sequence_id"]
    assert sequence_id == sample_sequence_id

    # Obtaining a pre mapping id from sample file
    sample_ga4gh_premapped_id = mapped_example["mapped_scores"][0]["pre_mapped"]["_id"]

    # Obtaining the same pre mapped ID from mapped variant, and its index
    for i in mapped_variant[urn][0]["pre_mapping"].keys():
        if mapped_variant[urn][0]["pre_mapping"][i]["_id"] == sample_ga4gh_premapped_id:
            index = i
            break

    # Obtaining post mapping ID from mapped variant
    computed_post = mapped_variant[urn][0]["mapped"][index]["_id"]

    # Obtaining post mapped ID from sample file
    sample_ga4gh_postmapped_id = mapped_example["mapped_scores"][0]["post_mapped"][
        "_id"
    ]

    assert computed_post == sample_ga4gh_postmapped_id
