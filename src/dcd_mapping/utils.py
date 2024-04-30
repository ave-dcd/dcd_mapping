"""Utility functions for dcd_mapping package"""
import json
from pathlib import Path
from typing import List, Optional

import hgvs.edit
import hgvs.location
import hgvs.parser
import hgvs.posedit
import hgvs.sequencevariant
from Bio.SeqUtils import seq3
from biocommons.seqrepo import SeqRepo
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import sha512t24u

from dcd_mapping.lookup import (
    get_chromosome_identifier,
    get_vrs_id_from_identifier,
)
from dcd_mapping.mavedb_data import get_raw_scoreset_metadata, get_scoreset_metadata
from dcd_mapping.schemas import (
    AlignmentResult,
    ComputedReferenceSequence,
    MappedOutput,
    MappedReferenceSequence,
    TargetSequenceType,
    TxSelectResult,
    VrsObject1_x,
)


def get_hgvs_string(allele: dict, dp: SeqRepo, ac: str) -> str:
    """Return an HGVS string for a given VRS allele

    :param allele: A post-mapped VRS allele
    :param dp: A SeqRepo instance
    :param acc: A RefSeq accession
    :return An HGVS string
    """
    stype = "p" if ac.startswith("NP") else "g"
    start = allele["location"]["interval"]["start"]["value"]
    end = allele["location"]["interval"]["end"]["value"]

    if start == end:
        ref = None
        aas = dp.get_sequence(ac, start - 1, start)
        aae = dp.get_sequence(ac, end, end + 1)
        end += 1
    else:
        ref = dp.get_sequence(ac, start, end)
        aas = dp.get_sequence(ac, start, start + 1)
        aae = dp.get_sequence(ac, end - 1, end)
        start += 1

    if stype == "p":
        ival = hgvs.location.Interval(
            start=hgvs.location.AAPosition(base=start, aa=aas),
            end=hgvs.location.AAPosition(base=end, aa=aae),
        )
    else:
        ival = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=start),
            end=hgvs.location.SimplePosition(base=end),
        )
    alt = allele["state"]["sequence"]

    edit = ""  # Set default
    if alt == ref:
        edit = "="
    if ref and (2 * ref == alt or len(ref) == 1 and set(ref) == set(alt)):
        edit = "dup"
    if alt == "":
        edit = "del"

    if edit != "dup" or edit != "del" or edit != "=":
        edit = (
            hgvs.edit.AARefAlt(ref=ref, alt=alt)
            if stype == "p"
            else hgvs.edit.NARefAlt(ref=ref, alt=alt)
        )

    if alt != ref:
        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
    else:
        posedit = f"{seq3(ref)}{start!s}=" if stype == "p" else f"{end!s}{ref}="

    var = str(hgvs.sequencevariant.SequenceVariant(ac=ac, type=stype, posedit=posedit))
    if var.endswith("delins"):
        var = var.replace("delins", "del")
    return var


def get_vod_premapped(allele: dict) -> dict:
    """Return a VariationDescriptor object given a VRS pre-mapped allele dict

    :param allele: A VRS allele dictionary
    :return A VariationDescriptor dictionary
    """
    return {
        "id": allele["id"],
        "type": "VariationDescriptor",
        "variation": {
            "id": allele["id"],
            "type": "Allele",
            "location": {
                "id": None,
                "type": "SequenceLocation",
                "sequence_id": allele["location"]["sequence_id"],
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "type": "Number",
                        "value": allele["location"]["interval"]["start"]["value"],
                    },
                    "end": {
                        "type": "Number",
                        "value": allele["location"]["interval"]["end"]["value"],
                    },
                },
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": allele["state"]["sequence"],
            },
        },
        "vrs_ref_allele_seq": allele["vrs_ref_allele_seq"],
    }


def get_vod_postmapped(allele: dict) -> dict:
    """Return a VariationDescriptor object given a VRS pre-mapped allele dict

    :param allele: A VRS allele dictionary
    :return A VariationDescriptor dictionary
    """
    return {
        "id": allele["id"],
        "type": "VariationDescriptor",
        "variation": {
            "id": allele["id"],
            "type": "Allele",
            "location": {
                "id": None,
                "type": "SequenceLocation",
                "sequence_id": allele["location"]["sequence_id"],
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "type": "Number",
                        "value": allele["location"]["interval"]["start"]["value"],
                    },
                    "end": {
                        "type": "Number",
                        "value": allele["location"]["interval"]["end"]["value"],
                    },
                },
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": allele["state"]["sequence"],
            },
        },
        "expressions": [
            {
                "type": "Expression",
                "syntax": "hgvs.p" if "p." in allele["hgvs"] else "hgvs.g",
                "value": allele["hgvs"],
                "syntax_version": None,
            }
        ],
        "vrs_ref_allele_seq": allele["vrs_ref_allele_seq"],
    }


def get_vod_haplotype(allele_list: List[dict]) -> dict:
    """Define VOD model for haplotype

    :param allele_list: A list of VRS allele dictionaries
    :return A VRS Haplotype-like structure
    """
    return {"type": "Haplotype", "members": allele_list}


def get_computed_reference_sequence(
    ss: str,
    layer: AnnotationLayer,
    tx_output: Optional[TxSelectResult] = None,
) -> ComputedReferenceSequence:
    """Report the computed reference sequence for a score set

    :param ss: A score set string
    :param layer: AnnotationLayer
    :param tx_output: Transcript data for a score set
    :return A ComputedReferenceSequence object
    """
    if layer == AnnotationLayer.PROTEIN:
        seq_id = f"ga4gh:SQ.{sha512t24u(tx_output.sequence.encode('ascii'))}"
        return ComputedReferenceSequence(
            sequence=tx_output.sequence,
            sequence_type=TargetSequenceType.PROTEIN,
            sequence_id=seq_id,
        )
    metadata = get_scoreset_metadata(ss)
    seq_id = f"ga4gh:SQ.{sha512t24u(metadata.target_sequence.encode('ascii'))}"
    return ComputedReferenceSequence(
        sequence=metadata.target_sequence,
        sequence_type=TargetSequenceType.DNA,
        sequence_id=seq_id,
    )


def get_mapped_reference_sequence(
    layer: AnnotationLayer,
    tx_output: Optional[TxSelectResult] = None,
    align_result: Optional[AlignmentResult] = None,
) -> MappedReferenceSequence:
    """Report the mapped reference sequence for a score set

    :param ss: A score set string
    :param layer: AnnotationLayer
    :param tx_output: Transcript data for a score set
    :return A MappedReferenceSequence object
    """
    if layer == AnnotationLayer.PROTEIN:
        return MappedReferenceSequence(
            sequence_type=TargetSequenceType.PROTEIN,
            sequence_id=get_vrs_id_from_identifier(tx_output.np),
            sequence_accessions=[tx_output.np],
        )
    seq_id = get_chromosome_identifier(align_result.chrom)
    return MappedReferenceSequence(
        sequence_type=TargetSequenceType.DNA,
        sequence_id=get_vrs_id_from_identifier(seq_id),
        sequence_accessions=[seq_id],
    )


def save_mapped_output_json(
    ss: str,
    mave_vrs_mappings: List[VrsObject1_x],
    align_result: AlignmentResult,
    tx_output: Optional[TxSelectResult] = None,
) -> None:
    """Save mapping output for a score set in a JSON file

    :param ss: Score set accession
    :param mave_vrs_mappings: A dictionary of VrsObject1_x objects
    :param align_result: Alignment information for a score set
    :tx_output: Transcript output for a score set
    :return None
    """
    curr = mave_vrs_mappings.variations
    layer = None
    for var in curr:
        if ss.startswith("urn:mavedb:00000097"):
            layer = AnnotationLayer.PROTEIN
            break
        if var.layer == AnnotationLayer.GENOMIC:
            layer = AnnotationLayer.GENOMIC
            break
    if not ss.startswith("urn:mavedb:00000097") and layer is None:
        layer = AnnotationLayer.PROTEIN

    mapped_ss_output = {}
    mapped_ss_output["metadata"] = get_raw_scoreset_metadata(ss)
    mapped_ss_output["computed_reference_sequence"] = get_computed_reference_sequence(
        ss=ss, layer=layer, tx_output=tx_output if tx_output else None
    ).model_dump()
    mapped_ss_output["mapped_reference_sequence"] = get_mapped_reference_sequence(
        tx_output=tx_output if tx_output else None,
        layer=layer,
        align_result=align_result,
    ).model_dump()

    mapped_scores = []
    for var in curr:
        if var and var.layer == layer:
            if "members" in var.pre_mapped_variants:
                pre_mapped_members = []
                post_mapped_members = []
                for sub_var in var.pre_mapped_variants["members"]:
                    pre_mapped_members.append(get_vod_premapped(sub_var))
                for sub_var in var.post_mapped_variants["members"]:
                    post_mapped_members.append(get_vod_postmapped(sub_var))
                mapped_scores.append(
                    MappedOutput(
                        pre_mapped=get_vod_haplotype(pre_mapped_members),
                        post_mapped=get_vod_haplotype(post_mapped_members),
                        mavedb_id=var.mavedb_id,
                        score=None if var.score == "NA" else float(var.score),
                    ).model_dump()
                )
            else:
                mapped_scores.append(
                    MappedOutput(
                        pre_mapped=get_vod_premapped(var.pre_mapped_variants),
                        post_mapped=get_vod_postmapped(var.post_mapped_variants),
                        mavedb_id=var.mavedb_id,
                        score=None if var.score == "NA" else float(var.score),
                    ).model_dump()
                )
    mapped_ss_output["mapped_scores"] = mapped_scores

    ss = ss.strip("urn:mavedb:")  # noqa: B005
    with (Path("analysis_files") / "mappings" / f"{ss}.json").open("w") as file:
        json.dump(mapped_ss_output, file, indent=4)
