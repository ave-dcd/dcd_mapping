import json
import pandas as pd

"""Function that specifies input format"""

def metadata_obtain(scoreset_json, scores_csv) -> dict:
    """
    Extract data from MaveDB scoresets and convert them into an input format that imitates
    the expected inputs of the package.


    Parameters
    ----------
        scoreset_json: json object
            Scoreset JSON object
        
        scores_csv

    Returns:
    ----------
        dat: dict
            Dictionary containing extracted metadata
        
        variant_data: dict
            Dictionary containing variants

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
    scores = vardat["score"].to_list()
    accessions = vardat["accession"].to_list()

    variant_data = {"hgvs_pro": varm, "hgvs_nt":ntlist, "scores":scores, "accessions":accessions}
    return dat, variant_data

#TODO: change other things according to these changes