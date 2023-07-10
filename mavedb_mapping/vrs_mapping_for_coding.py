# VRS Variant Mapping - Coding Scoresets
import io
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from gene.query import QueryHandler
from gene.database import create_db
from ga4gh.vrs.extras.translator import Translator
from Bio.Seq import Seq
from biocommons.seqrepo import SeqRepo
from transcript_selection_helper import *
import pandas as pd
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify, sha512t24u

sr = SeqRepo("/usr/local/share/seqrepo/latest", writeable=True)
dp = SeqRepoDataProxy(sr=sr)
tr = Translator(data_proxy=dp, normalize=False)
qh = QueryHandler(create_db())
vrs_mappings_dict = {}
scores_dict_coding = {}
mavedb_ids_coding = {}


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
        hgvs_string = ref + ":" + l + "." + varlist[i]
        allele = tr.translate_from(hgvs_string, "hgvs")

        if mapped == "pre":
            allele.location.sequence_id = "ga4gh:SQ." + sha512t24u(ts.encode("ascii"))
            if "dup" in hgvs_string:
                allele.state.sequence = 2 * str(
                    sr[str(allele.location.sequence_id)][
                        allele.location.start.value : allele.location.end.value
                    ]
                )

        else:
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
                        allele.location.interval.start.value = (
                            ranges[i][1] - diff - diff2
                        )
                        allele.location.interval.end.value = (
                            allele.location.interval.start.value + diff2
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
                    allele.location.interval.start.value : allele.location.interval.end.value
                ]
            )
        allele = normalize(allele, data_proxy=dp)
        allele._id = ga4gh_identify(allele)
        alleles.append(allele)

    if len(alleles) == 1:  # Not haplotype
        return alleles[0].as_dict()
    else:
        return models.Haplotype(members=alleles)


def process_protein_coding_data(
    varm, np, offset, tr, dp, ts, ranges, hits, scores, accessions
):
    var_ids_pre_map = []
    var_ids_post_map = []
    spro = []
    accpro = []

    for j in range(len(varm)):
        if (
            type(varm[j]) != str
            or len(varm[j]) == 3
            or varm[j] == "_wt"
            or varm[j] == "_sy"
        ):
            continue
        if varm[j].startswith("NP"):
            var_ids_pre_map.append(tr.translate_from(varm[j], "hgvs").as_dict())
            var_ids_post_map.append(tr.translate_from(varm[j], "hgvs").as_dict())
            spro.append(scores[j])
            accpro.append(accessions[j])
        else:
            try:
                if np.startswith("N") == True:
                    var_ids_pre_map.append(
                        get_haplotype_allele(
                            varm[j], np, 0, "p", tr, dp, ts, "pre", "", "", ""
                        )
                    )
                    var_ids_post_map.append(
                        get_haplotype_allele(
                            varm[j], np, offset, "p", tr, dp, ts, "post", "", "", ""
                        )
                    )
                    spro.append(scores[j])
                    accpro.append(accessions[j])
                else:
                    var_ids_pre_map.append(
                        get_haplotype_allele(
                            varm[j], np, 0, "p", tr, dp, ts, "pre", "", "", ""
                        )
                    )
                    var_ids_post_map.append(
                        get_haplotype_allele(
                            varm[j],
                            np,
                            offset,
                            "p",
                            tr,
                            dp,
                            ts,
                            "post",
                            ranges,
                            hits,
                            "",
                        )
                    )
                    spro.append(scores[j])
                    accpro.append(accessions[j])
            except:
                continue

    return var_ids_pre_map, var_ids_post_map, spro, accpro


def process_nt_data(ntlist, ref, ts, tr, dp, ranges, hits, scores, accessions, strand):
    var_ids_pre_map = []
    var_ids_post_map = []
    sn = []
    accn = []

    for j in range(len(ntlist)):
        if type(ntlist[j]) != str or ntlist[j] == "_wt" or ntlist[j] == "_sy":
            continue
        else:
            try:
                var_ids_pre_map.append(
                    get_haplotype_allele(
                        ntlist[j][2:],
                        ref,
                        0,
                        "g",
                        tr,
                        dp,
                        ts,
                        "pre",
                        ranges,
                        hits,
                        strand,
                    ).as_dict()
                )
                var_ids_post_map.append(
                    get_haplotype_allele(
                        ntlist[j][2:],
                        ref,
                        0,
                        "g",
                        tr,
                        dp,
                        ts,
                        "post",
                        ranges,
                        hits,
                        strand,
                    ).as_dict()
                )
                sn.append(scores[j])
                accn.append(accessions[j])
            except:
                continue

    return var_ids_pre_map, var_ids_post_map, sn, accn


def vrs_mapping_for_coding(dat, mappings_dict, mave_blat_dict, scores_csv):
    ranges = get_locs_list(mave_blat_dict["hits"])
    hits = get_hits_list(mave_blat_dict["hits"])
    ref = get_chr(dp, mave_blat_dict["chrom"])
    ts = dat["target_sequence"]
    if dat["target_type"] == "Protein coding" or dat["target_type"] == "protein_coding":
        item = mappings_dict
        vardat = pd.read_csv(io.StringIO(scores_csv.decode("utf-8")))
        scores = vardat["score"].to_list()
        accessions = vardat["accession"].to_list()
        varm = vardat["hgvs_pro"]

        mappings_list = []
        scores_list = []
        accessions_list = []

        np = mappings_dict["RefSeq_prot"]
        offset = mappings_dict["start"]
        ts = dat["target_sequence"]
        if len(set(str(ts))) > 4:
            stri = str(ts)
        else:
            ts = Seq(ts)
            ts = str(ts.translate(table=1)).replace("*", "")
        digest = "SQ." + sha512t24u(ts.encode("ascii"))
        alias_dict_list = [{"namespace": "ga4gh", "alias": digest}]
        sr.store(ts, nsaliases=alias_dict_list)

        var_ids_pre_map, var_ids_post_map, spro, accpro = process_protein_coding_data(
            varm, np, offset, tr, dp, ts, ranges, hits, scores, accessions
        )

        tempdat = pd.DataFrame(
            {"pre_mapping": var_ids_pre_map, "mapped": var_ids_post_map}
        )
        mappings_list.append(tempdat)
        scores_list.append(spro)
        accessions_list.append(accpro)

        if vardat["hgvs_nt"].isnull().values.all() == False and "97" not in dat["urn"]:
            ntlist = vardat["hgvs_nt"]
            varm = vardat["hgvs_pro"]
            strand = mave_blat_dict["strand"]
            var_ids_pre_map, var_ids_post_map, sn, accn = process_nt_data(
                ntlist, ref, ts, tr, dp, ranges, hits, scores, accessions, strand
            )

            tempdat = pd.DataFrame(
                {"pre_mapping": var_ids_pre_map, "mapped": var_ids_post_map}
            )
            mappings_list.append(tempdat)
            scores_list.append(sn)
            accessions_list.append(accn)

        vrs_mappings_dict[dat["urn"]] = mappings_list
        scores_dict_coding[dat["urn"]] = scores_list
        mavedb_ids_coding[dat["urn"]] = accessions_list

    return vrs_mappings_dict
