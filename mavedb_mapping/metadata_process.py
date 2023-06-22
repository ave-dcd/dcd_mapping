import requests
def metadata_obtain(urn):
    response = requests.get("https://api.mavedb.org/api/v1/score-sets/" + urn.replace(":","%3A"))
    json_parse = response.json()
    human_target_sequences = json_parse['targetGene']['wtSequence']['sequence']
    target_type = json_parse['targetGene']['wtSequence']['sequenceType']
    human_target = json_parse['targetGene']['name']
    uniprot = json_parse['targetGene']['externalIdentifiers'][0]['identifier']['identifier']
    target = json_parse['targetGene']['category']
    dat = {'urn': urn, 'target_sequence': human_target_sequences, 'target_sequence_type': target_type, 'target':human_target, 'uniprot_id':uniprot, 'target_type':target}
    return dat
 

urn = "urn:mavedb:00000041-a-1"
dat = metadata_obtain(urn)
print(dat)
