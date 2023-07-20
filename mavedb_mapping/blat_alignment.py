from Bio import SearchIO
import pandas as pd
import subprocess
from mavedb_mapping import qh
from mavedb_mapping import path_to_hg38_file

def get_gene_symb(dat):
    try:
        uniprot = dat["uniprot_id"]
        gsymb = qh.normalize(str(f"uniprot:{uniprot}")).gene_descriptor.label
    except:
        try:
            target = dat["target"].split(" ")[0]
            gsymb = qh.normalize(target).gene_descriptor.label
        except:
            return "NA"
    return gsymb


def get_gene_data(gsymb: str):
    if gsymb == "NA":
        return "NA"
    temp = qh.search(gsymb).source_matches
    source_dict = {}
    for i in range(len(temp)):
        source_dict[temp[i].source] = i
    return temp, source_dict


def get_hgnc_accession(temp, source_dict):
    accession = temp[source_dict["HGNC"]].records[0].concept_id
    return accession


def get_sequence_interval(records):
    for record in records:
        for location in record.locations:
            if location.interval.type == "SequenceInterval":
                start = location.interval.start.value
                end = location.interval.end.value
                loc_list = {"start": start, "end": end}
                return loc_list
    return None

def return_gene_data(return_chr: bool, temp, source_dict):
    """
    Parameters
    ----------
        return_chr :bool
           If True, returns chromosome information.

        dat: dict
            Dictionary containing data required for mapping.

        gsymb: str
            Gene symbol.


    Returns:
    --------
        str:
            If return_chr is True, returns the chromosome value as a string.
            If gene symbol cannot be extracted, returns 'NA'.
        OR
        dict:
            If gene symbol can be extracted and return_chr is False

    """
    if temp == "NA":
        return "NA"

    if "HGNC" in source_dict and return_chr == True:
        chrom = temp[source_dict["HGNC"]].records[0].locations[0].chr
        return chrom
    if source_dict.get("Ensembl") is not None and return_chr == False and len(temp[source_dict["Ensembl"]].records) != 0:
        loc_list = get_sequence_interval(temp[source_dict["Ensembl"]].records)
        if loc_list:
            return loc_list

    if source_dict.get("NCBI") is not None and return_chr == False and len(temp[source_dict["NCBI"]].records) != 0:
        loc_list = get_sequence_interval(temp[source_dict["NCBI"]].records)
        if loc_list:
            return loc_list

    return "NA"



def extract_blat_output(dat: dict):
    """
    Parameters
    ----------
        return_chr :bool
           If True, returns chromosome information.

        dat: dict
            Dictionary containing data required for mapping.


    Returns:
    --------
        str:
            If return_chr is True, returns the chromosome value as a string.
            If gene symbol cannot be extracted, returns 'NA'.
        OR
        dict:
            If gene symbol can be extracted and return_chr is False
    """
    blat_file = open("blat_query.fa", "w")
    blat_file.write(">" + dat["target"] + "\n")
    blat_file.write(dat["target_sequence"] + "\n")
    blat_file.close()
    if dat["target_sequence_type"] == "protein":
        command = (
            f"blat {path_to_hg38_file} -q=prot -t=dnax -minScore=20 blat_query.fa blat_out.psl"
        )
        process = subprocess.run(command, shell=True)
    else:
        command = f"blat {path_to_hg38_file} -minScore=20 blat_query.fa blat_out.psl"
        process = subprocess.run(command, shell=True)
    try:
        output = SearchIO.read("blat_out.psl", "blat-psl")
    except:
        try:
            command = (
                f"blat {path_to_hg38_file} -q=dnax -t=dnax -minScore=20 blat_query.fa blat_out.psl"
            )
            process = subprocess.run(command, shell=True)
            output = SearchIO.read("blat_out.psl", "blat-psl")
        except:
            return None
    return output


def get_query_and_hit_ranges(output, dat):
    """
    Extracts query and hit ranges from the BLAT output.

    Parameters
    ----------
        output:
            Output from the BLAT query.
        dat: dict
            Dictionary containing data required for mapping.

    Returns
    -------
        Tuple containing the chromosome, strand, coverage, identity, query ranges, hit ranges, and gene symbol.
    """
    hit_scores = list()
    hit_dict = {}
    use_chr = False
    chrom = strand = ""
    coverage = identity = None
    query_ranges = hit_ranges = list()

    gsymb = get_gene_symb(dat)
    temp, source_dict = get_gene_data(gsymb)
    accession = get_hgnc_accession(temp, source_dict)

    for c in range(len(output)):
        correct_chr = return_gene_data(True, temp, source_dict)
        if correct_chr == output[c].id.strip("chr"):
            use_chr = True
            break
        if (
            correct_chr == "NA"
        ):  # Take top scoring hit if target not found using gene normalizer
            hit_scores = list()
            for e in range(len(output[c])):
                hit_scores.append(output[c][e].score)
            hit_dict[c] = hit_scores
    if use_chr == False:
        for key in hit_dict:
            hit_dict[key] = max(hit_dict[key])
        hit = max(hit_dict, key=hit_dict.get)
    else:
        hit = c

    loc_dict = return_gene_data(False, temp, source_dict)

    hit_starts = list()
    for n in range(len(output[hit])):
        hit_starts.append(output[hit][n].hit_start)

    sub_scores = list()
    for n in range(len(output[hit])):
        sub_scores.append(output[hit][n].score)

    if loc_dict == "NA":
        hsp = output[hit][
            sub_scores.index(max(sub_scores))
        ]  # Take top score if no match found
    else:
        hsp = output[hit][
            hit_starts.index(min(hit_starts, key=lambda x: abs(x - loc_dict["start"])))
        ]

    for j in range(len(hsp)):
        test_file = open("blat_output_test.txt", "w")
        test_file.write(str(hsp[j]))
        test_file.close()

        query_string = ""
        hit_string = ""
        strand = hsp[0].query_strand
        coverage = 100 * (hsp.query_end - hsp.query_start) / output.seq_len
        # coverage = f"{hsp.query_end - hsp.query_start} / {output.seq_len}, {coverage}"
        identity = hsp.ident_pct

        test_file = open("blat_output_test.txt", "r")
        for k, line in enumerate(test_file):
            if k == 1:
                chrom = line.strip("\n")
            if k == 2:
                query_string = line.strip("\n")
            if k == 3:
                hit_string = line.strip("\n")
        test_file.close()

        chrom = chrom.split(" ")[9].strip("chr")
        query_string = query_string.split(" ")
        hit_string = hit_string.split(" ")
        query_ranges.append(query_string[2])
        hit_ranges.append(hit_string[4])

    return chrom, strand, coverage, identity, query_ranges, hit_ranges, gsymb, accession


def mave_to_blat(dat):
    """

    Performs BLAT Alignment on MaveDB scoreset data.

    Parameters
    ----------
        dat: dict
            Dictionary containing data from MaveDB scoresets.

    Returns
    -------
        mave_blat_dict: dict
            Dicitionary containing data after doing BLAT Alignment

    """
    output = extract_blat_output(dat)
    if output is not None:
        (
            chrom,
            strand,
            coverage,
            identity,
            query_ranges,
            hit_ranges,
            gsymb,
            accession,
        ) = get_query_and_hit_ranges(output, dat)
        qh_dat = {"query_ranges": query_ranges, "hit_ranges": hit_ranges}
        qh_dat = pd.DataFrame(data=qh_dat)
        mave_blat_dict = {
            "urn": dat["urn"],
            "chrom": chrom,
            "strand": strand,
            "target": dat["target"],
            "target_type": dat["target_type"],
            "uniprot": dat["uniprot_id"],
            "coverage": coverage,
            "identity": identity,
            "hits": qh_dat,
            "gsymb": gsymb,
            "accession": accession,
        }

    else:
        qh_dat = {"query_ranges": ["NA"], "hit_ranges": ["NA"]}
        qh_dat = pd.DataFrame(data=qh_dat)
        mave_blat_dict = {
            "urn": dat["urn"],
            "chrom": "NA",
            "strand": "NA",
            "target": "NA",
            "target_type": "NA",
            "uniprot": "NA",
            "coverage": "NA",
            "identity": "NA",
            "hits": qh_dat,
            "gsymb": "NA",
            "accession": "NA",
        }

    return mave_blat_dict
