import requests
from gene.query import QueryHandler
from gene.database import create_db
from ga4gh.core import ga4gh_identify, sha512t24u
from os import environ
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.normalize import normalize
from Bio.Seq import Seq
from ga4gh.vrs import models

sr = SeqRepo("/usr/local/share/seqrepo/latest", writeable=True)
environ["UTA_DB_URL"] = "postgresql://uta_admin:uta@localhost:5432/uta/uta_20210129"
qh = QueryHandler(create_db())


def get_start(string):
    return int(string.split(":")[0].strip("["))


def get_end(string):
    return int(string.split(":")[1].strip("]"))


def get_locs_list(hitsdat):
    locs_list = []
    for i in range(len(hitsdat.index)):
        start = get_start(hitsdat.at[i, "hit_ranges"])
        end = get_end(hitsdat.at[i, "hit_ranges"])
        locs_list.append([start, end])
    return locs_list


def get_hits_list(hitsdat):
    hits_list = []
    for i in range(len(hitsdat.index)):
        start = get_start(hitsdat.at[i, "query_ranges"])
        end = get_end(hitsdat.at[i, "query_ranges"])
        hits_list.append([start, end])
    return hits_list


def get_query_hits(dat):
    query_list = []
    hits_list = []
    for i in range(len(dat.index)):
        query_start = get_start(dat.at[i, "query_ranges"])
        query_end = get_end(dat.at[i, "query_ranges"])
        query_list.append([query_start, query_end])
        hit_start = get_start(dat.at[i, "hit_ranges"])
        hit_end = get_end(dat.at[i, "hit_ranges"])
        hits_list.append([hit_start, hit_end])
        return query_list, hits_list


def get_ga4gh(dp, ref):
    aliases = dp.get_metadata(ref)["aliases"]
    f = filter(lambda x: "ga4gh" in x, aliases)
    return "ga4gh:" + list(f)[0].split(":")[1]


def get_chr(dp, chrom):
    aliases = dp.get_metadata("GRCh38:" + chrom)["aliases"]
    f = filter(lambda x: "refseq" in x, aliases)
    return list(f)[0].split(":")[1]


def modify_hgvs(var, ref, off, hp):
    if len(var) == 3 or var == "_wt" or var == "_sy" or "[" in var:
        return var
    var = ref + ":" + var
    var = hp.parse_hgvs_variant(var)
    var.posedit.pos.start.base = var.posedit.pos.start.base + off
    var.posedit.pos.end.base = var.posedit.pos.end.base + off
    return str(var)


def blat_check(mave_blat_dict, dat):
    if mave_blat_dict["uniprot"] == None:
        test = dat["target"].split(" ")
        for j in range(len(test)):
            try:
                out = qh.normalize(test[j]).gene_descriptor
                gene_dat = [out.label, out.extensions[2].value["chr"]]
                if mave_blat_dict["chrom"] != gene_dat[1]:
                    return False
                else:
                    return True
            except:
                continue


def get_haplotype_allele(var, ref, offset, l, tr, dp, ts, mapped, ranges, hits, strand):
    var = var.lstrip(f"{l}.")

    if "[" in var:
        var = var[1:][:-1]
        varlist = var.split(";")
        varlist = list(set(varlist))
    else:
        varlist = list()
        varlist.append(var)

    locs = {}
    alleles = []

    for i in range(len(varlist)):
        try:
            hgvs_string = ref + ":" + l + "." + varlist[i]
            allele = tr.translate_from(hgvs_string, "hgvs")

            if mapped == "pre":
                allele.location.sequence_id = "ga4gh:SQ." + sha512t24u(
                    ts.encode("ascii")
                )
                if "dup" in hgvs_string:
                    allele.state.sequence = 2 * str(
                        sr[str(allele.location.sequence_id)][
                            allele.location.start.value : allele.location.end.value
                        ]
                    )

            else:
                if l != "g":
                    allele.location.start.value = allele.location.start.value + offset
                    allele.location.end.value = allele.location.end.value + offset
                    if "dup" in hgvs_string:
                        allele.state.sequence = 2 * str(
                            sr[str(allele.location.sequence_id)][
                                allele.location.start.value : allele.location.end.value
                            ]
                        )

                else:
                    start = allele.location.start.value
                    if len(hits) == 1 and strand == 1:
                        i = 0
                        diff = start - hits[i][0]
                        diff2 = allele.location.end.value - start
                        allele.location.start.value = ranges[i][0] + diff
                        allele.location.end.value = allele.location.start.value + diff2
                    else:
                        for i in range(len(hits)):
                            if start >= hits[i][0] and start < hits[i][1]:
                                break
                        diff = start - hits[i][0]
                        diff2 = allele.location.end.value - start
                        if strand == 1:  # positive orientation
                            allele.location.start.value = ranges[i][0] + diff
                            allele.location.end.value = (
                                allele.location.start.value + diff2
                            )
                            if "dup" in hgvs_string:
                                allele.state.sequence = 2 * str(
                                    sr[str(allele.location.sequence_id)][
                                        allele.location.start.value : allele.location.end.value
                                    ]
                                )
                        else:
                            allele.location.start.value = ranges[i][1] - diff - diff2
                            allele.location.end.value = (
                                allele.location.start.value + diff2
                            )
                            if "dup" in hgvs_string:
                                allele.state.sequence = 2 * str(
                                    sr[str(allele.location.sequence_id)][
                                        allele.location.start.value : allele.location.end.value
                                    ]
                                )
                            allele.state.sequence = str(
                                Seq(str(allele.state.sequence)).reverse_complement()
                            )

            if allele.state.sequence == "N" and l != "p":
                allele.state.sequence = str(
                    sr[str(allele.location.sequence_id)][
                        allele.location.start.value : allele.location.end.value
                    ]
                )
            allele = normalize(allele, data_proxy=dp)
            allele.id = ga4gh_identify(allele)
            alleles.append(allele)
        except:
            vrstext = {"definition": ref + ":" + l + "." + varlist[i], "type": "Text"}
            return vrstext

    if len(alleles) == 1:  # Not haplotype
        return alleles[0]
    else:
        return models.Haplotype(members=alleles)


def get_clingen_id(hgvs):
    url = "https://reg.genome.network/allele?hgvs=" + hgvs
    page = requests.get(url).json()
    page = page["@id"]
    try:
        return page.split("/")[4]
    except:
        return "NA"
