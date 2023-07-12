from Bio import SearchIO
from gene.query import QueryHandler
from Bio import SearchIO
import pandas as pd
import subprocess
from gene.database import create_db


class BLATALignment:
    def __init__(self, dat):
        self.qh = QueryHandler(create_db())
        self.dat = dat

    def get_gene_symb(self):
        try:
            uniprot = self.dat["uniprot_id"]
            gsymb = self.qh.normalize(str(f"uniprot:{uniprot}")).gene_descriptor.label
        except:
            try:
                target = self.dat["target"].split(" ")[0]
                gsymb = self.qh.normalize(target).gene_descriptor.label
            except:
                return "NA"
        return gsymb

    def get_gene_data(self, return_chr: bool, gsymb: str):
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
        if gsymb == "NA":
            return "NA"
        temp = self.qh.search(gsymb).source_matches
        source_dict = {}
        for i in range(len(temp)):
            source_dict[temp[i].source] = i

        hgnc_accesion = temp[source_dict["HGNC"]].records[0].concept_id

        if "HGNC" in source_dict and return_chr == True:
            chrom = temp[source_dict["HGNC"]].records[0].locations[0].chr
            return chrom, hgnc_accesion

        if (
            "Ensembl" in source_dict
            and return_chr == False
            and len(temp[source_dict["Ensembl"]].records) != 0
        ):
            for j in range(len(temp[source_dict["Ensembl"]].records)):
                for k in range(len(temp[source_dict["Ensembl"]].records[j].locations)):
                    if (
                        temp[source_dict["Ensembl"]]
                        .records[j]
                        .locations[k]
                        .interval.type
                        == "SequenceInterval"
                    ):  # Multiple records per source
                        start = (
                            temp[source_dict["Ensembl"]]
                            .records[j]
                            .locations[k]
                            .interval.start.value
                        )
                        end = (
                            temp[source_dict["Ensembl"]]
                            .records[j]
                            .locations[k]
                            .interval.end.value
                        )
                        loc_list = {}
                        loc_list["start"] = start
                        loc_list["end"] = end
                        return loc_list, hgnc_accesion
        if (
            "NCBI" in source_dict
            and return_chr == False
            and len(temp[source_dict["NCBI"]].records) != 0
        ):
            for j in range(len(temp[source_dict["NCBI"]].records)):
                for k in range(len(temp[source_dict["NCBI"]].records[j].locations)):
                    if (
                        temp[source_dict["NCBI"]].records[j].locations[k].interval.type
                        == "SequenceInterval"
                    ):
                        start = (
                            temp[source_dict["NCBI"]]
                            .records[j]
                            .locations[k]
                            .interval.start.value
                        )
                        end = (
                            temp[source_dict["NCBI"]]
                            .records[j]
                            .locations[k]
                            .interval.end.value
                        )
                        loc_list = {}
                        loc_list["start"] = start
                        loc_list["end"] = end
                        return loc_list, hgnc_accesion
        return "NA"

    def extract_blat_output(self):
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
        blat_file.write(">" + self.dat["target"] + "\n")
        blat_file.write(self.dat["target_sequence"] + "\n")
        blat_file.close()
        if self.dat["target_sequence_type"] == "protein":
            command = (
                "blat hg38.2bit -q=prot -t=dnax -minScore=20 blat_query.fa blat_out.psl"
            )
            process = subprocess.run(command, shell=True)
        else:
            command = "blat hg38.2bit -minScore=20 blat_query.fa blat_out.psl"
            process = subprocess.run(command, shell=True)
        try:
            output = SearchIO.read("blat_out.psl", "blat-psl")
        except:
            try:
                command = "blat hg38.2bit -q=dnax -t=dnax -minScore=20 blat_query.fa blat_out.psl"
                process = subprocess.run(command, shell=True)
                output = SearchIO.read("blat_out.psl", "blat-psl")
            except:
                return None
        return output

    def get_query_and_hit_ranges(self, output):
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
        gsymb = self.get_gene_symb()

        for c in range(len(output)):
            correct_chr, hgnc_accesion = self.get_gene_data(
                return_chr=True, gsymb=gsymb
            )
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

        loc_dict, hgnc_accesion = self.get_gene_data(False, gsymb)

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
                hit_starts.index(
                    min(hit_starts, key=lambda x: abs(x - loc_dict["start"]))
                )
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

        return (
            chrom,
            strand,
            coverage,
            identity,
            query_ranges,
            hit_ranges,
            gsymb,
            hgnc_accesion,
        )

    def check_non_human(mave_blat_dict, min_percentage=80):
        # for dna min % = 95
        # for prot min % = 80
        # as per BLAT website: "BLAT on DNA is designed to quickly find sequences of 95% and grent or shorter sequence alignments. BLAT on proteins finds sequences of 80% and greater similarity of length 20 amino acids or more.
        """
        Checks if a sample is human or non-human based on the Mave-BLAT dictionary.

        Parameters
        ----------

            mave_blat_dict: dict
                Dicitionary containing data after doing BLAT Alignment

            min_percent: int
                Minimum percentage coverage to consider a sample as human.

        Returns
        -------
            str: "human" if the sample is human, "Non human" otherwise.
        """
        cov = mave_blat_dict["coverage"]
        if cov == "NA":
            return "Non human"
        else:
            if cov < min_percentage:
                return "Non human"
            else:
                return "human"

    def mave_to_blat(self):
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
        output = self.extract_blat_output()
        if output is not None:
            (
                chrom,
                strand,
                coverage,
                identity,
                query_ranges,
                hit_ranges,
                gsymb,
                hgnc_accesion,
            ) = self.get_query_and_hit_ranges(output)
            qh_dat = {"query_ranges": query_ranges, "hit_ranges": hit_ranges}
            qh_dat = pd.DataFrame(data=qh_dat)
            mave_blat_dict = {
                "urn": self.dat["urn"],
                "chrom": chrom,
                "strand": strand,
                "target": self.dat["target"],
                "target_type": self.dat["target_type"],
                "uniprot": self.dat["uniprot_id"],
                "coverage": coverage,
                "identity": identity,
                "hits": qh_dat,
                "gsymb": gsymb,
                "hgnc_accesion": hgnc_accesion,
            }

        else:
            qh_dat = {"query_ranges": ["NA"], "hit_ranges": ["NA"]}
            qh_dat = pd.DataFrame(data=qh_dat)
            mave_blat_dict = {
                "urn": self.dat["urn"],
                "chrom": "NA",
                "strand": "NA",
                "target": "NA",
                "target_type": "NA",
                "uniprot": "NA",
                "coverage": "NA",
                "identity": "NA",
                "hits": qh_dat,
                "gsymb": "NA",
            }

        return mave_blat_dict
