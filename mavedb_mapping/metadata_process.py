import json
def metadata_obtain(scoreset_json) -> dict:
    '''
    Parameters
    ----------
        scoreset_json: json object
            Scoreset JSON object
    
    Returns:
    ----------
        dat: dict
            Metadata extracted from the scoreset JSON.

    '''
    scoreset = json.load(scoreset_json)
    urn = scoreset["urn"]

    #Sequence
    human_target_sequences = scoreset['targetGene']['wtSequence']['sequence']

    #Sequence type
    target_type = scoreset['targetGene']['wtSequence']['sequenceType']

    #Name of gene/protein
    human_target = scoreset['targetGene']['name']

    #Uniprot ID
    if scoreset['targetGene']['externalIdentifiers']!=[]:
        uniprot = scoreset['targetGene']['externalIdentifiers'][0]['identifier']['identifier']
    else:
        uniprot = None

    #Organism name
    organism = scoreset['targetGene']['referenceMaps'][0]['genome']['organismName']

    #Target type
    target = scoreset['targetGene']['category']


    dat = {'urn': urn, 'target_sequence': human_target_sequences, 'target_sequence_type': target_type, 'target':human_target,  'uniprot_id':uniprot, 'target_type':target, 'organism': organism}
    return dat
 
