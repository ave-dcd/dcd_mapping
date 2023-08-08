import json


def metadata_obtain(scoreset_json) -> dict:
    """
    Extracts relevant metadata from a scoreset JSON object
    Parameters
    ----------
        scoreset_json: json object
            Scoreset JSON object

    Returns:
    ----------
        dat: dict
            Dictionary containing extracted metadata
    """
    scoreset = json.load(scoreset_json)
    urn = scoreset["urn"]

    # Sequence
    human_target_sequences = scoreset["targetGene"]["wtSequence"]["sequence"]

    # Sequence type
    target_type = scoreset["targetGene"]["wtSequence"]["sequenceType"]

    # Uniprot ID
    if scoreset["targetGene"]["externalIdentifiers"] != []:
        uniprot = scoreset["targetGene"]["externalIdentifiers"][0]["identifier"][
            "identifier"
        ]
    else:
        uniprot = None

    # Target type
    target = scoreset["targetGene"]["category"]

    dat = {
        "urn": urn,
        "target_sequence": human_target_sequences,
        "target_sequence_type": target_type,
        "uniprot_id": uniprot,
        "target_type": target,
    }
    return dat
