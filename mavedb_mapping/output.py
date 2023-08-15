import json


def output(dat, mappings, vrsmaps, blat_dict):
    output_dict = {}
    if mappings:
        output_dict["sequence_metadata"] = coding(dat, mappings, blat_dict)
        output_dict["mapping"] = mapped_variants(vrsmaps[0])
        return output_dict
    else:
        # TODO:for non coding
        pass


def coding(dat, mappings, blat_dict):
    seq_dict = {
        "Seq": dat["target_sequence"],
        "RefSeq_prot": mappings["RefSeq_prot"],
        "RefSeq_nuc": mappings["RefSeq_nuc"],
        "gsymb": mappings["gsymb"],
        "chrom": blat_dict["chrom"],
    }
    return seq_dict


def mapped_variants(vrs):
    premapped = vrs["pre_mapping"]
    mapped = vrs["mapped"]
    mapp = list()
    print(type(mapped))
    for i in range(len(mapped)):
        c = maps(premapped[i], mapped[i])
        mapp.append(c)
    return mapp


def maps(pre, post, relation="SO:is_homologous_to"):
    d = {}
    d["pre_mapped"] = pre
    d["mapped"] = post
    d["relation"] = relation
    return d
