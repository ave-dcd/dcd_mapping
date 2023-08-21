# Find start location in provided target sequence when start position is not first position of sequence
import Bio
import re
from Bio.Seq import Seq
from Bio.SeqUtils import seq1


# maybe put in get haplotype
def validation_helper(protstring):
    protstring = protstring[1:][:-1]
    vs = protstring.split(";")
    return vs


def is_true(i, locs, ts, aa_dict):
    for j in range(len(locs)):
        if ts[i + locs[j] - locs[0]] != aa_dict[locs[j]]:
            return False
    return True


def if_start_not_first(dat: dict, vardat):
    if dat["target_type"] == "Protein coding" and dat["target_sequence_type"] == "dna":
        oseq = Seq(dat["target_sequence"])
        ts = str(oseq.translate(table=1))

        protlist = vardat["hgvs_pro"]

        if type(protlist[0]) != str or protlist[0].startswith("NP"):
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

            seq = ""
            for key in aa_dict:
                seq = seq + aa_dict[key]

            for i in range(len(ts)):
                if is_true(i, locs, ts, aa_dict):
                    if i + 1 == min(aa_dict.keys()) or i + 2 == min(aa_dict.keys()):
                        offset = 0
                    else:
                        offset = i
                    break
            return offset
