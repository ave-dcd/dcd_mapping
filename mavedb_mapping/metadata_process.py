import json
def metadata_obtain(scoreset_json):
    scoreset = json.load(scoreset_json)
    urn = scoreset["urn"]
    human_target_sequences = scoreset['targetGene']['wtSequence']['sequence']
    target_type = scoreset['targetGene']['wtSequence']['sequenceType']
    if scoreset['targetGene']['externalIdentifiers']!=[]:
        uniprot = scoreset['targetGene']['externalIdentifiers'][0]['identifier']['identifier']
    else:
        uniprot = None
    target = scoreset['targetGene']['category']
    dat = {'urn': urn, 'target_sequence': human_target_sequences, 'target_sequence_type': target_type,'uniprot_id':uniprot, 'target_type':target}
    return dat
 
