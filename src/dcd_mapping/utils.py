"""HGVS utility function for dcd_mapping package"""

import hgvs.edit
import hgvs.location
import hgvs.parser
import hgvs.posedit
import hgvs.sequencevariant
from Bio.SeqUtils import seq3
from biocommons.seqrepo import SeqRepo


def get_hgvs_string(allele: dict, dp: SeqRepo, ac: str) -> str:
    """Return an HGVS string for a given VRS allele
    :param allele: A post-mapped VRS allele
    :param dp: A SeqRepo instance
    :param acc: A RefSeq accession
    :return An HGVS string
    """
    if ac.startswith('NP'):
        stype = 'p'
    else:
        stype = 'g'
    start = allele["location"]["interval"]["start"]["value"]
    end = allele["location"]["interval"]["end"]["value"]

    if start == end:
        ref = None
        aas = dp.get_sequence(ac, start-1, start)
        aae = dp.get_sequence(ac, end, end + 1)
        end += 1
    else:
        ref = dp.get_sequence(ac, start, end)
        aas = dp.get_sequence(ac, start, start + 1)
        aae = dp.get_sequence(ac, end - 1, end)
        start += 1

    if stype == 'p':
        ival = hgvs.location.Interval(start=hgvs.location.AAPosition
                                      (base=start, aa = aas),
                                      end=hgvs.location.AAPosition(base=end, aa = aae))
    else:
        ival = hgvs.location.Interval(start=hgvs.location.SimplePosition(
            base=start), end=hgvs.location.SimplePosition(base=end))
    alt = allele["state"]["sequence"]

    edit = '' # Set default
    if alt == ref:
        edit = '='
    if ref:
        if 2*ref == alt or len(ref) == 1 and set(ref) == set(alt):
            edit = 'dup'
    if alt == '':
        edit = 'del'

    if edit != 'dup' or edit != 'del' or edit != '=':
        if stype == 'p':
            edit = hgvs.edit.AARefAlt(ref=ref, alt=alt)
        else:
            edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)

    if alt != ref:
        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
    else:
        if stype == 'p':
            posedit = seq3(ref) + str(start) + '='
        else:
            posedit = str(end) + ref + '='

    var = str(hgvs.sequencevariant.SequenceVariant(ac=ac,type=stype,posedit=posedit))
    if var.endswith('delins'):
        var = var.replace('delins', 'del')
    return var
