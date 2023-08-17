from ga4gh.vrs.extras.translator import Translator
from Bio.Seq import Seq
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs.normalize import normalize
from mavedb_mapping import sr, dp

tr = Translator(data_proxy=dp, normalize=False)


def return_normalized_allele(allele, l):
    if allele.state.sequence == "N" and l != "p":
        allele.state.sequence = str(
            sr[str(allele.location.sequence_id)][
                allele.location.interval.start.value : allele.location.interval.end.value
            ]
        )
    allele = normalize(allele, data_proxy=dp)
    allele._id = ga4gh_identify(allele)
    return allele


def pre_mapping_allele_normalize(ref, var, l, ts):
    hgvs_string = ref + ":" + l + "." + var
    allele = tr.translate_from(hgvs_string, "hgvs")
    allele.location.sequence_id = "ga4gh:SQ." + sha512t24u(ts.encode("ascii"))
    if "dup" in hgvs_string:
        allele.state.sequence = 2 * str(
            sr[str(allele.location.sequence_id)][
                allele.location.start.value : allele.location.end.value
            ]
        )
    allele = return_normalized_allele(allele, l)
    return allele


def pre_mapping(ref, var, l, ts):
    var = var.lstrip(f"{l}.")
    if "[" in var:
        var = var[1:][:-1]
        varlist = var.split(";")
        varlist = list(set(varlist))
    else:
        allele = pre_mapping_allele_normalize(ref, var, l, ts)
        return allele.as_dict()

    alleles = []
    for i in varlist:
        allele = pre_mapping_allele_normalize(ref, i, l, ts)
        alleles.append(allele.as_dict())
    return models.Haplotype(members=alleles)


def post_mapping(ref, var, l, ranges, hits, offset, strand):
    var = var.lstrip(f"{l}.")
    if "[" in var:
        var = var[1:][:-1]
        varlist = var.split(";")
        varlist = list(set(varlist))
    else:
        allele = post_mapping_allele_return(ref, var, l, ranges, hits, offset, strand)
        return allele
    alleles = []
    for i in varlist:
        allele = post_mapping_allele_return(ref, i, l, ranges, hits, offset, strand)
        alleles.append(allele)
    return models.Haplotype(members=alleles)


def post_mapping_allele_return(ref, var, l, ranges, hits, offset, strand):
    hgvs_string = ref + ":" + l + "." + var
    allele = tr.translate_from(hgvs_string, "hgvs")
    if l != "g":
        allele.location.interval.start.value += offset
        allele.location.interval.end.value += offset
        if "dup" in hgvs_string:
            allele.state.sequence = 2 * str(
                sr[str(allele.location.sequence_id)][
                    allele.location.interval.start.value : allele.location.interval.end.value
                ]
            )
    else:
        start = allele.location.interval.start.value
        if len(hits) == 1 and strand == 1:
            i = 0
            diff = start - hits[i][0]
            diff2 = allele.location.interval.end.value - start
            allele.location.interval.start.value = ranges[i][0] + diff
            allele.location.interval.end.value = (
                allele.location.interval.start.value + diff2
            )
        else:
            for i in range(len(hits)):
                if start >= hits[i][0] and start < hits[i][1]:
                    break
            diff = start - hits[i][0]
            diff2 = allele.location.interval.end.value - start
            if strand == 1:  # positive orientation
                allele.location.interval.start.value = ranges[i][0] + diff
                allele.location.interval.end.value = (
                    allele.location.interval.start.value + diff2
                )
                if "dup" in hgvs_string:
                    allele.state.sequence = 2 * str(
                        sr[str(allele.location.sequence_id)][
                            allele.location.interval.start.value : allele.location.interval.end.value
                        ]
                    )
            else:
                allele.location.interval.start.value = ranges[i][1] - diff - diff2
                allele.location.interval.end.value = (
                    allele.location.interval.start.value + diff2
                )
                if "dup" in hgvs_string:
                    allele.state.sequence = 2 * str(
                        sr[str(allele.location.sequence_id)][
                            allele.location.interval.start.value : allele.location.interval.end.value
                        ]
                    )
                    allele.state.sequence = str(
                        Seq(str(allele.state.sequence)).reverse_complement()
                    )
    return allele
