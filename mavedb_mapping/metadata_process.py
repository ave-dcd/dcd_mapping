import json
def metadata_obtain(scoreset_json):
    scoreset = json.load(scoreset_json)
    urn = scoreset["urn"]
    human_target_sequences = scoreset['targetGene']['wtSequence']['sequence']
    target_type = scoreset['targetGene']['wtSequence']['sequenceType']
    human_target = scoreset['targetGene']['name']
    human_assembly = 0000
    if scoreset['targetGene']['externalIdentifiers']!=[]:
        uniprot = scoreset['targetGene']['externalIdentifiers'][0]['identifier']['identifier']
    else:
        uniprot = None
    organism = scoreset['targetGene']['referenceMaps'][0]['genome']['organismName']
    target = scoreset['targetGene']['category']
    dat = {'urn': urn, 'target_sequence': human_target_sequences, 'target_sequence_type': target_type, 'target':human_target,  'uniprot_id':uniprot, 'target_type':target, 'organism': organism}
    return dat
 
