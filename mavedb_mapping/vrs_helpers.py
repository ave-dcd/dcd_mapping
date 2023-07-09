from ga4gh.vrs.normalize import normalize
from Bio.Seq import Seq
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify, sha512t24u
sr = SeqRepo("/usr/local/share/seqrepo/latest", writeable = True)



vrs_mappings_dict = {}
scores_dict_coding = {}
mavedb_ids_coding = {}


def get_haplotype_allele(var, ref, offset, l, tr, dp, ts, mapped, ranges, hits, strand):
    var = var.lstrip(f'{l}.')

    if '[' in var:
        var = var[1:][:-1]
        varlist = var.split(';')
        varlist = list(set(varlist))
    else:
        varlist = list()
        varlist.append(var)
  
    alleles = []

    for varitem in varlist:
        hgvs_string = f'{ref}:{l}.{varitem}'
        allele = tr.translate_from(hgvs_string, 'hgvs')
        try:
            if mapped == 'pre':
                allele.location.sequence_id = f'ga4gh:SQ.{sha512t24u(ts.encode("ascii"))}'
                if 'dup' in hgvs_string:
                    allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])
                
            else:
                if l != 'g':
                    allele.location.interval.start.value += offset
                    allele.location.interval.end.value += offset
                    if 'dup' in hgvs_string:
                        allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.interval.start.value:allele.location.interval.end.value])
                        
                else:
                    start = allele.location.interval.start.value
                    if len(hits) == 1 and strand == 1:
                        diff = start - hits[0][0]
                        diff2 = allele.location.interval.end.value - start
                        allele.location.interval.start.value = ranges[i][0] + diff
                        allele.location.interval.end.value = allele.location.interval.start.value + diff2
                    else:
                        for i in range(len(hits)):
                            if start >= hits[i][0] and start < hits[i][1]:
                                break
                        diff = start - hits[i][0]
                        diff2 = allele.location.interval.end.value - start
                        if strand == 1: # positive orientation
                            allele.location.interval.start.value = ranges[i][0] + diff
                            allele.location.interval.end.value = allele.location.interval.start.value + diff2
                            if 'dup' in hgvs_string:
                                allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.interval.start.value:allele.location.interval.end.value])
                        else: 
                            allele.location.interval.start.value = ranges[i][1] - diff - diff2
                            allele.location.interval.end.value = allele.location.interval.start.value + diff2
                            if 'dup' in hgvs_string:
                                allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])
                            allele.state.sequence = str(Seq(str(allele.state.sequence)).reverse_complement())
            
            if allele.state.sequence == 'N' and l != 'p':
                allele.state.sequence = str(sr[str(allele.location.sequence_id)][allele.location.interval.start.value:allele.location.interval.end.value])
            allele = normalize(allele, data_proxy = dp)    
            allele._id = ga4gh_identify(allele)
            alleles.append(allele)
        except:
            vrstext = {'definition':ref + ':'+ l +'.' + varlist[i], 'type': 'Text'}
            return vrstext
       
    if len(alleles) == 1: # Not haplotype
        return alleles[0].as_dict()
    else:
        return models.Haplotype(members = alleles)
    

def add_custom_digest_to_seqrepo(ts):
    if len(set(str(ts))) > 4:
        stri = str(ts)  
    else:
        ts = Seq(ts)
        ts = str(ts.translate(table=1)).replace('*', '')
    digest = 'SQ.' + sha512t24u(ts.encode('ascii'))
    alias_dict_list = [{'namespace': 'ga4gh', 'alias': digest}]
    sr.store(ts, nsaliases = alias_dict_list)
