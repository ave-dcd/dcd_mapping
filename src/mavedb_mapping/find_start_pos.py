# Find start location in provided target sequence when start position is not first position of sequence
import requests
import Bio
import re
from Bio.Seq import Seq
import pandas as pd
import io
from Bio.SeqUtils import seq1
from Bio.Seq import Seq

offset_within_ts = {}


def validation_helper(protstring):
    protstring = protstring[1:][:-1]
    vs = protstring.split(";")
    return vs


def if_start_not_first(dat):
    if dat["target_type"] == "Protein coding" and dat["target_sequence_type"] == "dna":
        urn = dat["urn"]
        if (
            urn == "urn:mavedb:00000053-a-1" or urn == "urn:mavedb:00000053-a-2"
        ):  # target sequence missing codon
            return None
        oseq = Seq(dat["target_sequence"])
        ts = str(oseq.translate(table=1))

        string = (
            "https://api.mavedb.org/api/v1/score-sets/urn%3Amavedb%3A"
            + dat["urn"][11::]
            + "/scores"
        )
        # TODO: remove api call
        origdat = requests.get(string).content
        score_dat = pd.read_csv(io.StringIO(origdat.decode("utf-8")))
        protlist = score_dat["hgvs_pro"].to_list()
        if type(score_dat.at[0, "hgvs_pro"]) != str or score_dat.at[
            0, "hgvs_pro"
        ].startswith("NP"):
            return None
        protlist = [x.lstrip("p.") for x in protlist]

        aa_dict = {}

        for k in range(len(protlist)):
            if protlist[k] == "_sy" or protlist[k] == "_wt":
                continue
            else:
                if ";" in protlist[k]:
                    vs = validation_helper(protlist[k])
                    for l in range(len(vs)):
                        aa = vs[l][:3]
                        if (
                            aa == "="
                            or vs[l][-3:]
                            not in Bio.SeqUtils.IUPACData.protein_letters_3to1.keys()
                        ):
                            continue
                        if "=" in vs[l]:
                            loc = vs[l][3:][:-1]
                        else:
                            loc = vs[l][3:][:-3]
                        if loc not in aa_dict:
                            loc = re.sub("[^0-9]", "", loc)
                            aa_dict[loc] = seq1(aa)

                else:
                    if "_" in protlist[k]:
                        continue
                    aa = protlist[k][:3]
                    if (
                        aa == "="
                        or protlist[k][-3:]
                        not in Bio.SeqUtils.IUPACData.protein_letters_3to1.keys()
                    ):
                        continue
                    if "=" in protlist[k]:
                        loc = protlist[k][3:][:-1]
                    else:
                        loc = protlist[k][3:][:-3]
                    if loc not in aa_dict:
                        loc = re.sub("[^0-9]", "", loc)
                        aa_dict[loc] = seq1(aa)

        aa_dict.pop("", None)

        err_locs = []
        for m in range(len(ts)):
            if str(m) in list(aa_dict.keys()):
                if aa_dict[str(m)] != ts[int(m) - 1]:  # Str vs dict offset
                    err_locs.append(m)

        if len(err_locs) > 1:
            aa_dict = {int(k): v for k, v in aa_dict.items()}
            aa_dict = sorted(aa_dict.items())
            aa_dict = dict(aa_dict)
            locs = list(aa_dict.keys())[0:5]
            p0, p1, p2, p3, p4 = locs[0], locs[1], locs[2], locs[3], locs[4]
            offset = locs[0]

            seq = ""
            for key in aa_dict:
                seq = seq + aa_dict[key]

            for i in range(len(ts)):
                if (
                    ts[i] == aa_dict[p0]
                    and ts[i + p1 - p0] == aa_dict[p1]
                    and ts[i + p2 - p0] == aa_dict[p2]
                    and ts[i + p3 - p0] == aa_dict[p3]
                    and ts[i + p4 - p0] == aa_dict[p4]
                ):
                    if i + 1 == min(aa_dict.keys()) or i + 2 == min(aa_dict.keys()):
                        offset_within_ts[urn] = 0
                    else:
                        offset_within_ts[urn] = i
                    break
        return offset_within_ts
