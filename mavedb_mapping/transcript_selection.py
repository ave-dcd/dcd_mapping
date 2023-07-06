import requests
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from gene.query import QueryHandler
import nest_asyncio
import asyncio
from cool_seq_tool.data_sources.uta_database import UTADatabase
from cool_seq_tool.data_sources.mane_transcript_mappings import MANETranscriptMappings
from os import environ
from Bio.Seq import Seq
from biocommons.seqrepo import SeqRepo
from bs4 import BeautifulSoup
from transcript_selection_helper import *

sr = SeqRepo("/usr/local/share/seqrepo/latest", writeable=True)
environ["UTA_DB_URL"] = "postgresql://uta_admin:uta@localhost:5432/uta/uta_20210129"
utadb = UTADatabase(db_pwd="uta")
mane = MANETranscriptMappings()
qh = QueryHandler(create_db())
dp = SeqRepoDataProxy(sr=sr)


nest_asyncio.apply()


async def mapq(locs, chrom, gsymb):
    transcript_lists = []
    for i in range(len(locs)):
        testquery = f"""select *
                            from uta_20210129.tx_exon_aln_v
                            where hgnc = '{gsymb}'
                            and {locs[i][0]} between alt_start_i and alt_end_i
                            or {locs[i][1]} between alt_start_i and alt_end_i
                            and alt_ac = '{chrom}'"""

        out = await utadb.execute_query(testquery)
        tl = []
        for j in range(len(out)):
            if out[j]["tx_ac"].startswith("NR_") == False:
                tl.append(out[j]["tx_ac"])
        if tl != []:
            transcript_lists.append(tl)
    return transcript_lists


def get_gsymb(dat):
    try:
        uniprot = dat["uniprot_id"]
        gsymb = qh.normalize(str(f"uniprot:{uniprot}")).gene_descriptor.label
    except:
        temp = dat["target"].split(" ")
        if temp[0] == "JAK":
            temp[0] = "JAK1"
        gsymb = qh.normalize(temp[0]).gene_descriptor.label
    return gsymb


def using_uniprot(dat: dict):
    try:  # Look for transcripts using uniprot id
        url = "https://www.uniprot.org/uniprot/" + str(dat["uniprot_id"]) + ".xml"
        page = requests.get(url)
        page = BeautifulSoup(page.text)
        page = page.find_all("sequence")
        up = page[1].get_text()
        stri = str(dat["target_sequence"])
        if up.find(stri) != -1:
            full_match = True
        else:
            full_match = False
        start = up.find(stri[:10])
        return start, full_match
    except:
        # print(dat['urn'], 'no transcripts found')
        return "NA", "NA"


def get_status(mane_trans: list):
    if len(mane_trans) == 1:
        np = mane_trans[0]["RefSeq_prot"]
        nm = mane_trans[0]["RefSeq_nuc"]
        status = "MANE Select"
    else:
        if mane_trans[0]["MANE_status"] == "MANE Select":
            np = mane_trans[0]["RefSeq_prot"]
            nm = mane_trans[0]["RefSeq_nuc"]
            status = "MANE Select"
        else:
            np = mane_trans[1]["RefSeq_prot"]
            nm = mane_trans[1]["RefSeq_nuc"]
            status = "MANE Plus Clinical"
    return status, np, nm


def from_mane_trans(dat, mane_trans):
    oseq = dat["target_sequence"]
    status, np, nm = get_status(mane_trans)
    if len(set(str(oseq))) > 4:
        stri = str(oseq)
    else:
        oseq = Seq(oseq)
        stri = str(oseq.translate(table=1)).replace("*", "")

    if str(sr[np]).find(stri) != -1:
        full_match = True
    else:
        full_match = False
    start = str(sr[np]).find(stri[:10])
    mappings_list = [np, start, dat["urn"], full_match, nm, status]
    return np, start, full_match, nm, status


async def np(nm):
    testquery = (
        f"SELECT pro_ac FROM uta_20210129.associated_accessions WHERE tx_ac = '{nm}'"
    )
    out = await utadb.execute_query(testquery)
    try:
        return out[0]["pro_ac"]
    except:
        return out


def no_mane_trans(isect, dat):
    trans_lens = []
    for i in range(len(isect)):
        trans_lens.append(len(str(sr[isect[i]])))
        loc = trans_lens.index(max(trans_lens))
        nm = isect[loc]
        np = asyncio.run(np(nm))

    if np != []:
        oseq = dat["target_sequence"]

        if len(set(str(oseq))) > 4:
            stri = str(oseq)
        else:
            oseq = Seq(oseq)
            stri = str(oseq.translate(table=1)).replace("*", "")

        if str(sr[np]).find(stri) != -1:
            full_match = True
        else:
            full_match = False
        start = str(sr[np]).find(stri[:10])
        status = "Longest Compatible"
        return np, start, full_match, nm, status


def main(mave_blat_dict, dat):
    if dat["target_type"] == "Protein coding" or dat["target_type"] == "protein_coding":
        if mave_blat_dict["chrom"] == "NA":
            return "NA"
        locs = get_locs_list(mave_blat_dict["hits"])
        chrom = get_chr(dp, mave_blat_dict["chrom"])
        gsymb = get_gsymb(dat)
        ts = asyncio.run(mapq(locs, chrom, gsymb))
        try:
            isect = list(set.intersection(*map(set, ts)))
        except:
            start, full_match = using_uniprot(dat)
            np = nm = status = "NA"
            mappings_dict = {
                "urn": dat["urn"],
                "uniprot_id": dat["uniprot_id"],
                "start": start,
                "full_match": full_match,
                "RefSeq_prot": "NA",
                "RefSeq_nuc": "NA",
                "status": "NA",
                "gsymb": gsymb,
            }
            return mappings_dict

        mane_trans = mane.get_mane_from_transcripts(isect)
        if mane_trans != []:
            np, start, full_match, nm, status = from_mane_trans(dat, mane_trans)
        else:
            np, start, full_match, nm, status = no_mane_trans(isect, dat)

        mappings_dict = {
            "urn": dat["urn"],
            "uniprot_id": dat["uniprot_id"],
            "start": start,
            "full_match": full_match,
            "RefSeq_prot": np,
            "RefSeq_nuc": nm,
            "status": status,
            "gsymb": gsymb,
            "uniprot_id": dat["uniprot_id"],
        }
        return mappings_dict
