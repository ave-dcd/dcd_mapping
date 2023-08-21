from ga4gh.vrs.extras.translator import Translator
from Bio.Seq import Seq
from mavedb_mapping.transcript_selection_helper import HelperFunctionsForBLATOutput
import pandas as pd
from ga4gh.core import sha512t24u
from mavedb_mapping import sr, dp
from mavedb_mapping.get_allele import pre_mapping, post_mapping

tr = Translator(data_proxy=dp, normalize=False)

vrs_mappings_dict = {}
scores_dict_coding = {}
mavedb_ids_coding = {}


def process_protein_coding_data(
    varm: list,
    np: str,
    offset: int,
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
            HGVS protein variant list from MaveDB scores

        np:str
            RefSeq Protein Identifier.

        offset: int
            Offset value.

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
        tempdat: dict
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
                    var_ids_pre_map.append(pre_mapping(np, varm[j], "p", ts))
                    var_ids_post_map.append(
                        post_mapping(np, varm[j], "p", ranges, hits, offset, None)
                    )
                    spro.append(scores[j])
                    accpro.append(accessions[j])
                else:
                    var_ids_pre_map.append(pre_mapping(np, varm[j], "p", ts))

                    var_ids_post_map.append(
                        post_mapping(np, varm[j], "p", ranges, hits, offset, None)
                    )
                    spro.append(scores[j])
                    accpro.append(accessions[j])
            except:
                continue
    tempdat = {"pre_mapping": var_ids_pre_map, "mapped": var_ids_post_map}
    return tempdat, spro, accpro


def process_nt_data(ntlist, ref, ts, ranges, hits, scores, accessions, strand):
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

        ranges: list
            List of ranges from BLAT

        hits: list
            List of hits from BLAT

        scores: list
            List of scores from MaveDB

        accessions: list
            List of accessions from MaveDB scores


    Returns
    -------
        tempdat: dict
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
                var_ids_pre_map.append(pre_mapping(ref, ntlist[j][2:], "g", ts))
                var_ids_post_map.append(
                    post_mapping(ref, ntlist[j][2:], "g", ranges, hits, 0, strand)
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
    varm = variant_data["hgvs_pro"]  # TODO: change variable name l
    ntlist = variant_data["hgvs_nt"]

    np = mappings_dict["RefSeq_prot"]
    offset = mappings_dict["start"]

    ts = dat["target_sequence"]
    ts = Seq(ts)
    ts = str(ts.translate(table=1)).replace("*", "")
    digest = "SQ." + sha512t24u(ts.encode("ascii"))
    alias_dict_list = [{"namespace": "ga4gh", "alias": digest}]
    sr.store(ts, nsaliases=alias_dict_list)
    print(2)
    mappings_list = list()
    scores_list = list()
    accessions_list = list()

    tempdat, spro, accpro = process_protein_coding_data(
        varm, np, offset, ts, ranges, hits, scores, accessions
    )

    mappings_list.append(tempdat)
    scores_list.append(spro)
    accessions_list.append(accpro)

    if ntlist.isnull().values.all() == False:
        strand = mave_blat_dict["strand"]
        tempdat, sn, accn = process_nt_data(
            ntlist, ref, ts, ranges, hits, scores, accessions, strand
        )

        mappings_list.append(tempdat)
        scores_list.append(sn)
        accessions_list.append(accn)

    return mappings_list


def check_for_transcripts(mappings_dict):
    if mappings_dict["status"] == "NA":
        raise Exception("No transcripts found")


def vrs_mapping_for_non_coding(dat, mave_blat_dict, ref, ranges, hits, variant_data):
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
        ntlist, ref, ts, ranges, hits, scores, accessions, strand
    )

    mappings_list.append(tempdat)
    scores_list.append(sn)
    accessions_list.append(accn)

    return mappings_list


def vrs_mapping(
    dat: dict, mappings_dict: dict, mave_blat_dict: dict, variant_data: dict
) -> list:
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

        variant_data: dict
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
            dat, mave_blat_dict, ref, ranges, hits, variant_data
        )

    return mapping
