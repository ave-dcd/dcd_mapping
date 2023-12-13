"""Produces the output with all data from the mapping"""


def output(dat: dict, mappings: dict, vrsmaps: dict, blat_dict: dict) -> dict:
    """
    Produces final output in the form of a dictionary

    Parameters
    ----------
        dat: dict
            Dictionary containing data required for mapping.

        mappings: dict
            Dictionary after transcript selection.

        vrsmaps: dict
            Dictionary containing premapped and mapped variants

        blat_dict: dict
            Dicitionary containing data after doing BLAT Alignment
    """
    output_dict = {}
    if mappings:
        output_dict["sequence"] = protein_coding(dat, mappings, blat_dict)
        output_dict["mapping"] = mapped_variants(vrsmaps[0])
        if len(vrsmaps > 1):
            output_dict["mapping"] = mapped_variants(vrsmaps[1])

        return output_dict
    else:
        output_dict["sequence"] = non_coding(dat, blat_dict)
        output_dict["mapping"] = mapped_variants(vrsmaps[0])


def non_coding(dat, blat_dict, relation="SO:is_homologous_to"):
    seq_dict = {
        "target": dat["target_sequence"],
        "chrom": blat_dict["chrom"],
        "target_type": dat["target_type"],
        "relation": relation,
    }
    return seq_dict


def protein_coding(dat, mappings, blat_dict, relation="SO:is_homologous_to"):
    seq_dict = {
        "target": dat["target_sequence"],
        "RefSeq_prot": mappings["RefSeq_prot"],
        "RefSeq_nuc": mappings["RefSeq_nuc"],
        "gsymb": mappings["gsymb"],
        "chrom": blat_dict["chrom"],
        "target_type": dat["target_type"],
        "relation": relation,
    }
    return seq_dict


def mapped_variants(vrs):
    premapped = vrs["pre_mapping"]
    mapped = vrs["mapped"]
    mapp = list()
    for i in range(len(premapped)):
        c = maps(premapped[i], mapped[i])
        mapp.append(c)
    return mapp


def maps(pre, post):
    d = {}
    d["pre_mapped"] = pre
    d["mapped"] = post

    return d
