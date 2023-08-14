import json
import pandas as pd


def metadata_obtain(scoreset_json, scores_csv) -> dict:
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

    # Sequence
    human_target_sequences = scoreset["targetGene"]["wtSequence"]["sequence"]

    # Sequence type
    target_type = scoreset["targetGene"]["wtSequence"]["sequenceType"]

    # Target type
    target = scoreset["targetGene"]["category"]

    dat = {
        "target_sequence": human_target_sequences,
        "target_sequence_type": target_type,
        "target_type": target,
    }
    
    vardat = pd.read_csv(scores_csv)

    varm = vardat["hgvs_pro"]
    ntlist = vardat["hgvs_nt"]

    variant_data = {"hgvs_pro": varm, "hgvs_nt":ntlist}
    return dat, variant_data

#TODO: change other things according to these changes