class HelperFunctionsForBLATOutput:
    def __init__(self, mave_blat_dict):
        self.mave_blat_dict = mave_blat_dict
        self.hitsdat = mave_blat_dict["hits"]

    @staticmethod
    def get_start(string):
        return int(string.split(":")[0].strip("["))

    @staticmethod
    def get_end(string):
        return int(string.split(":")[1].strip("]"))

    def get_locs_list(self):
        locs_list = []
        for i in range(len(self.hitsdat.index)):
            start = self.get_start(self.hitsdat.at[i, "hit_ranges"])
            end = self.get_end(self.hitsdat.at[i, "hit_ranges"])
            locs_list.append([start, end])
        return locs_list

    def get_hits_list(self):
        hits_list = []
        for i in range(len(self.hitsdat.index)):
            start = self.get_start(self.hitsdat.at[i, "query_ranges"])
            end = self.get_end(self.hitsdat.at[i, "query_ranges"])
            hits_list.append([start, end])
        return hits_list

    def get_chr(self, dp):
        aliases = dp.get_metadata("GRCh38:" + self.mave_blat_dict["chrom"])["aliases"]
        f = filter(lambda x: "refseq" in x, aliases)
        return list(f)[0].split(":")[1]

    def is_human(self, min_percentage=80):
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
            bool: True if the sample is human, False otherwise.
        """
        cov = self.mave_blat_dict["coverage"]
        if cov == "NA":
            return False
        else:
            if cov < min_percentage:
                return False
            else:
                return True
