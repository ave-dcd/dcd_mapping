from ga4gh.vrs.extras.translator import Translator
from Bio.Seq import Seq
from mavedb_mapping.transcript_selection_helper import HelperFunctionsForBLATOutput
import pandas as pd
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs.normalize import normalize
from mavedb_mapping import sr, dp

tr = Translator(data_proxy=dp, normalize=False)

vrs_mappings_dict = {}
scores_dict_coding = {}
mavedb_ids_coding = {}


def get_haplotype_allele(
    var: str,
    ref: str,
    offset: int,
    l: str,
    tr,
    ts: str,
    mapped: str,
    ranges,
    hits,
    strand,
):
    """
    Retrieves allele, and normalizes it

    Parameters
    ----------
        var: str
            Variant String.

        ref:str
            RefSeq Identifier.

        offset: int
            Offset Value.

        l:str

        tr: GA4GH Translator.

        dp: Data proxy.

        ts:str
            Target Sequence.

        mapped: str
            Mapping state: premapping/ postmapping.

        ranges: list
            List of ranges

        hits: list

        strand: int


    Returns:
    --------
        allele: dict
            Normalized Allele

    """

    var = var.lstrip(f"{l}.")

    if "[" in var:
        var = var[1:][:-1]
        varlist = var.split(";")
        varlist = list(set(varlist))
    else:
        varlist = list()
        varlist.append(var)
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
    varm: list,
    np: str,
    offset: int,
    tr,
    ts,
    ranges,
    hits,
    scores: list,
    accessions: list,
):
    """
    Uses HGVS protein data to map alleles

    Parameters
    ----------
        varm: list
            HGVS protein list from MaveDB scores

        np:str
            RefSeq Protein Identifier.

        offset: int
            Offset value.

        tr: GA4GH Translator

        dp: Data proxy

        ts: str
            Target Sequence

        ranges: list
            List of ranges

        hits: list
            List of hits

        scores: list
            List of scores from MaveDB

        accessions: list
            List of accessions from MaveDB scores


    Returns
    -------
        tempdat: DataFrame
            Contains premapped and postmapped alleles

        spro: list
            Scores List

        accpro: list
            Accessions List

    """
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
                            varm[j], np, 0, "p", tr, ts, "pre", "", "", ""
                        )
                    )
                    var_ids_post_map.append(
                        get_haplotype_allele(
                            varm[j], np, offset, "p", tr, ts, "post", "", "", ""
                        )
                    )
                    spro.append(scores[j])
                    accpro.append(accessions[j])
                else:
                    var_ids_pre_map.append(
                        get_haplotype_allele(
                            varm[j], np, 0, "p", tr, ts, "pre", "", "", ""
                        )
                    )
                    var_ids_post_map.append(
                        get_haplotype_allele(
                            varm[j],
                            np,
                            offset,
                            "p",
                            tr,
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
    tempdat = pd.DataFrame({"pre_mapping": var_ids_pre_map, "mapped": var_ids_post_map})
    return tempdat, spro, accpro


def process_nt_data(ntlist, ref, ts, tr, ranges, hits, scores, accessions, strand):
    """
    Uses HGVS nucleotide data to map alleles

    Parameters
    ----------
        ntlist: list
            HGVS nucleotide list from MaveDB scores

        ref: str
            reference sequence accession number (NCBI)

        ts: str
            Target Sequence

        tr: GA4GH Translator

        dp: Data proxy

        ranges: list
            List of ranges

        hits: list
            List of hits

        scores: list
            List of scores from MaveDB

        accessions: list
            List of accessions from MaveDB scores


    Returns
    -------
        tempdat: DataFrame
            Contains premapped and postmapped alleles

        sn: list
            Scores List

        accn: list
            Accessions List

    """
    var_ids_pre_map = []
    var_ids_post_map = []
    sn = []
    accn = []

    for j in range(len(ntlist)):
        if ntlist[j] == "_wt" or ntlist[j] == "_sy":
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
                        ts,
                        "pre",
                        ranges,
                        hits,
                        strand,
                    )
                )
                var_ids_post_map.append(
                    get_haplotype_allele(
                        ntlist[j][2:],
                        ref,
                        0,
                        "g",
                        tr,
                        ts,
                        "post",
                        ranges,
                        hits,
                        strand,
                    )
                )
                sn.append(scores[j])
                accn.append(accessions[j])
            except:
                continue
    tempdat = pd.DataFrame({"pre_mapping": var_ids_pre_map, "mapped": var_ids_post_map})
    return tempdat, sn, accn


def vrs_mapping_for_protein_coding(
    dat, mappings_dict, mave_blat_dict, ref, ranges, hits, variant_data
):
    """
    Perform VRS mapping for protein coding scoresets.

    Parameters
    ----------
        dat: dict
            Dictionary containing data required for mapping.

        mappings_dict: dict
            Dictionary after transcript selection.

        mave_blat_dict: dict
            Dicitionary containing data after doing BLAT Alignment

        ref: str
            reference sequence accession number (NCBI)

        ranges: list
            List of ranges.

        hits: list
            List of hits.

        scores_csv: csv
            Scores from MaveDB.

    Returns
    -------
        vrs_mappings_dict: dict
            VRS mappings dictionary.
    """
    scores = variant_data["scores"]
    accessions = variant_data["accessions"]
    varm = variant_data["hgvs_pro"] #TODO: change variable name l
    ntlist = variant_data["hgvs_nt"]

    np = mappings_dict["RefSeq_prot"]
    offset = mappings_dict["start"]

    ts = dat["target_sequence"]
    ts = Seq(ts)
    ts = str(ts.translate(table=1)).replace("*", "")
    digest = "SQ." + sha512t24u(ts.encode("ascii"))
    alias_dict_list = [{"namespace": "ga4gh", "alias": digest}]
    sr.store(ts, nsaliases=alias_dict_list)

    mappings_list = list()
    scores_list = list()
    accessions_list = list()

    tempdat, spro, accpro = process_protein_coding_data(
        varm, np, offset, tr, ts, ranges, hits, scores, accessions
    )

    mappings_list.append(tempdat)
    scores_list.append(spro)
    accessions_list.append(accpro)

    if ntlist.isnull().values.all() == False and "97" not in dat["urn"]:
        strand = mave_blat_dict["strand"]
        tempdat, sn, accn = process_nt_data(
            ntlist, ref, ts, tr, ranges, hits, scores, accessions, strand
        )

        mappings_list.append(tempdat)
        scores_list.append(sn)
        accessions_list.append(accn)

    vrs_mappings_dict[dat["urn"]] = mappings_list
    scores_dict_coding[dat["urn"]] = scores_list
    mavedb_ids_coding[dat["urn"]] = accessions_list

    return vrs_mappings_dict


def check_for_transcripts(mappings_dict):
    if mappings_dict["status"] == "NA":
        raise Exception("No transcripts found")


def vrs_mapping_for_non_coding(
    dat, mappings_dict, mave_blat_dict, ref, ranges, hits, variant_data
):
    """
    Perform VRS mapping for protein coding scoresets.

    Parameters
    ----------
        dat: dict
            Dictionary containing data required for mapping.

        mappings_dict: dict
            Dictionary after transcript selection.

        mave_blat_dict: dict
            Dicitionary containing data after doing BLAT Alignment

        ref: str
            reference sequence accession number (NCBI)

        ranges: list
            List of ranges.

        hits: list
            List of hits.

        scores_csv: csv
            Scores from MaveDB.

    Returns
    -------
        vrs_mappings_dict: dict
            VRS mappings dictionary.
    """

    mappings_list = list()
    scores_list = list()
    accessions_list = list()
    scores = variant_data["scores"]
    accessions = variant_data["accessions"]
    ntlist = variant_data["hgvs_nt"]
    ts = Seq(dat["target_sequence"])
    ts = str(ts.translate(table=1)).replace("*", "")
    digest = "SQ." + sha512t24u(ts.encode("ascii"))
    alias_dict_list = [{"namespace": "ga4gh", "alias": digest}]
    sr.store(ts, nsaliases=alias_dict_list)

    strand = mave_blat_dict["strand"]
    tempdat, sn, accn = process_nt_data(
        ntlist, ref, ts, tr, ranges, hits, scores, accessions, strand
    )

    mappings_list.append(tempdat)
    scores_list.append(sn)
    accessions_list.append(accn)

    vrs_mappings_dict[dat["urn"]] = mappings_list
    scores_dict_coding[dat["urn"]] = scores_list
    mavedb_ids_coding[dat["urn"]] = accessions_list

    return vrs_mappings_dict


def vrs_mapping(dat, mappings_dict, mave_blat_dict, variant_data):
    """
    Perform VRS mapping for protein coding scoresets.

    Parameters
    ----------
        dat: dict
            Dictionary containing data required for mapping.

        mappings_dict: dict
            Dictionary after transcript selection.

        mave_blat_dict: dict
            Dicitionary containing data after doing BLAT Alignment

        scores_csv: csv
            Scores from MaveDB.

    Returns
    -------
        mapping: dict
            VRS mappings dictionary.
    """
    helper = HelperFunctionsForBLATOutput(mave_blat_dict)
    ranges = helper.get_locs_list()
    hits = helper.get_hits_list()
    ref = helper.get_chr()

    if dat["target_type"] == "Protein coding" and dat["target_sequence_type"] == "dna":
        check_for_transcripts(mappings_dict)
        mapping = vrs_mapping_for_protein_coding(
            dat, mappings_dict, mave_blat_dict, ref, ranges, hits, variant_data
        )
    else:
        mapping = vrs_mapping_for_non_coding(
            dat, mappings_dict, mave_blat_dict, ref, ranges, hits, variant_data
        )

    return mapping
