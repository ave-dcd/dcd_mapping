# prot and non prot
import pytest
import json
from mavedb_mapping import data_file_path

def get_sample_mapping_data(urn):
    """Function to obtain data from sample mapping"""

    # Sample Mapping file
    mapped_file = open(f"{data_file_path}{urn[11:]}.json")
    mapped_example = json.load(mapped_file)
    print(mapped_example.keys())

    # Sequence id from sample mapping file
    sample_sequence_id = mapped_example["mapped_scores"][0]["pre_mapped"]["location"][
        "sequence_id"
    ]

    # Obtaining a pre mapping id from sample file
    sample_ga4gh_premapped_id = mapped_example["mapped_scores"][0]["pre_mapped"]["_id"]

    # Obtaining post mapped ID from sample file
    sample_ga4gh_postmapped_id = mapped_example["mapped_scores"][0]["post_mapped"][
        "_id"
    ]

    return sample_sequence_id, sample_ga4gh_premapped_id, sample_ga4gh_postmapped_id


def test_for_vrs_mapping(full_mapping):
    """Test function to compare results from sample mappings obtained from the notebooks"""

    mapped_variant, urn = full_mapping

    (
        sample_sequence_id,
        sample_ga4gh_premapped_id,
        sample_ga4gh_postmapped_id,
    ) = get_sample_mapping_data(urn)

    # Sequence ID from mapped variant
    sequence_id = mapped_variant[urn][0]["pre_mapping"][0]["location"]["sequence_id"]
    assert sequence_id == sample_sequence_id

    # Obtaining the same pre mapped ID from mapped variant, and its index
    for i in mapped_variant[urn][0]["pre_mapping"].keys():
        if mapped_variant[urn][0]["pre_mapping"][i]["_id"] == sample_ga4gh_premapped_id:
            index = i
            break

    # Obtaining post mapping ID from mapped variant
    computed_post = mapped_variant[urn][0]["mapped"][index]["_id"]

    assert computed_post == sample_ga4gh_postmapped_id
