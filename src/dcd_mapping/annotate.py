"""Annotate MaveDB score set metadata with mapped scores."""
import json
import logging
from pathlib import Path

import hgvs.edit
import hgvs.location
import hgvs.parser
import hgvs.posedit
import hgvs.sequencevariant
from Bio.SeqUtils import seq3
from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import sha512t24u

from dcd_mapping.lookup import (
    get_chromosome_identifier,
    get_chromosome_identifier_from_vrs_id,
    get_seqrepo,
    get_vrs_id_from_identifier,
)
from dcd_mapping.mavedb_data import get_raw_scoreset_metadata, get_scoreset_metadata
from dcd_mapping.resource_utils import LOCAL_STORE_PATH
from dcd_mapping.schemas import (
    AlignmentResult,
    ComputedReferenceSequence,
    MappedOutput,
    MappedReferenceSequence,
    ScoresetMetadata,
    TargetSequenceType,
    TxSelectResult,
    VrsMapping1_3,
)

_logger = logging.getLogger(__name__)


def get_hgvs_string(allele: dict, dp: SeqRepoAccess, ac: str) -> str:
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


def get_vod_haplotype(allele_list: list[dict]) -> dict:
    """Define VOD model for haplotype

    :param allele_list: A list of VRS allele dictionaries
    :return A VRS Haplotype-like structure
    """
    return {"type": "Haplotype", "members": allele_list}


def get_computed_reference_sequence(
    ss: str,
    layer: AnnotationLayer,
    tx_output: TxSelectResult | None = None,
) -> ComputedReferenceSequence:
    """Report the computed reference sequence for a score set

    :param ss: A score set string
    :param layer: AnnotationLayer
    :param tx_output: Transcript data for a score set
    :return A ComputedReferenceSequence object
    """
    if layer == AnnotationLayer.PROTEIN:
        if tx_output is None:
            raise ValueError
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
    tx_output: TxSelectResult | None = None,
    align_result: AlignmentResult | None = None,
) -> MappedReferenceSequence:
    """Report the mapped reference sequence for a score set

    :param ss: A score set string
    :param layer: AnnotationLayer
    :param tx_output: Transcript data for a score set
    :return A MappedReferenceSequence object
    """
    if layer == AnnotationLayer.PROTEIN and tx_output is not None:
        vrs_id = get_vrs_id_from_identifier(tx_output.np)
        if vrs_id is None:
            raise ValueError
        return MappedReferenceSequence(
            sequence_type=TargetSequenceType.PROTEIN,
            sequence_id=vrs_id,
            sequence_accessions=[tx_output.np],
        )
    seq_id = get_chromosome_identifier(align_result.chrom)
    vrs_id = get_vrs_id_from_identifier(seq_id)
    if vrs_id is None:
        raise ValueError
    return MappedReferenceSequence(
        sequence_type=TargetSequenceType.DNA,
        sequence_id=vrs_id,
        sequence_accessions=[seq_id],
    )


def _set_layer(ss: str, mappings: list[VrsMapping1_3]) -> AnnotationLayer:
    if ss.startswith("urn:mavedb:00000097"):
        return AnnotationLayer.PROTEIN
    for var in mappings:
        if var.layer == AnnotationLayer.GENOMIC:
            return AnnotationLayer.GENOMIC
    return AnnotationLayer.PROTEIN


def _format_score_mapping(var: VrsMapping1_3, layer: AnnotationLayer) -> dict | None:
    if var and var.layer == layer:
        if "members" in var.pre_mapped_variants:
            pre_mapped_members = []
            post_mapped_members = []
            for sub_var in var.pre_mapped_variants["members"]:
                pre_mapped_members.append(get_vod_premapped(sub_var))
            for sub_var in var.post_mapped_variants["members"]:
                post_mapped_members.append(get_vod_postmapped(sub_var))
            return MappedOutput(
                pre_mapped=get_vod_haplotype(pre_mapped_members),
                post_mapped=get_vod_haplotype(post_mapped_members),
                mavedb_id=var.mavedb_id,
                score=None if var.score == "NA" else float(var.score),
            ).model_dump()
        return MappedOutput(
            pre_mapped=get_vod_premapped(var.pre_mapped_variants),
            post_mapped=get_vod_postmapped(var.post_mapped_variants),
            mavedb_id=var.mavedb_id,
            score=None if var.score == "NA" else float(var.score),
        ).model_dump()
    return None


def save_mapped_output_json(
    urn: str,
    mappings: list[VrsMapping1_3],
    align_result: AlignmentResult,
    tx_output: TxSelectResult | None = None,
    output_path: Path | None = None,
) -> None:
    """Save mapping output for a score set in a JSON file

    :param urn: Score set accession
    :param mave_vrs_mappings: A dictionary of VrsObject1_x objects
    :param align_result: Alignment information for a score set
    :param tx_output: Transcript output for a score set
    :param output_path:
    """
    layer = _set_layer(urn, mappings)

    mapped_ss_output = {
        "metadata": get_raw_scoreset_metadata(urn),
        "computed_reference_sequence": get_computed_reference_sequence(
            ss=urn, layer=layer, tx_output=tx_output
        ).model_dump(),
        "mapped_reference_sequence": get_mapped_reference_sequence(
            tx_output=tx_output,
            layer=layer,
            align_result=align_result,
        ).model_dump(),
    }

    mapped_scores = []
    for var in mappings:
        formatted_score_mapping = _format_score_mapping(var, layer)
        mapped_scores.append(formatted_score_mapping)
    mapped_ss_output["mapped_scores"] = mapped_scores

    urn = urn.removeprefix("urn:mavedb:")
    if not output_path:
        output_path = LOCAL_STORE_PATH / f"{urn}_mapping.json"

    with output_path.open("w") as file:
        json.dump(mapped_ss_output, file, indent=4)


def _format_start_end(ss: str, start: int, end: int) -> list[int]:
    """Format start and end coordinates for vrs_ref_allele_seq for known edge cases

    :param ss: score set
    :param start: start coordinate
    :param end: end coordinate
    :return A list of start and end coordinates
    """
    if ss.startswith("urn:mavedb:00000060-a-1"):
        _logger.warning(
            "urn:mavedb:00000060-a-1 reports the entire human reference sequence as the target sequence, but the start and end positions need to be manually offset by 289"
        )
        return [start + 289, end + 289]
    if ss.startswith("urn:mavedb:00000060-a-2"):
        _logger.warning(
            "urn:mavedb:00000060-a-2 reports the entire human reference sequence as the target sequence, but the start and end positions need to be manually offset by 331"
        )
        return [start + 331, end + 331]
    return [start, end]


def annotate(
    tx_select_results: TxSelectResult | None,
    vrs_results: list[VrsMapping1_3],
    metadata: ScoresetMetadata,
) -> list[VrsMapping1_3]:
    """Given a list of mappings, add additional contextual data:

    1. ``vrs_ref_allele_seq``: The sequence between the start and end positions
        indicated in the variant
    2. ``hgvs``: An HGVS string describing the variant (only included for post-mapped
        variants)

    :param tx_select_results: transcript selection if available
    :param vrs_results: in-progress variant mappings
    :param metadata: MaveDB scoreset metadata
    :return: annotated mappings objects
    """
    sr = get_seqrepo()
    for var in vrs_results:
        if not var:
            continue
        variant_list = var.pre_mapped_variants
        if "members" in variant_list:
            for sub_var in variant_list["members"]:
                start_end = _format_start_end(
                    metadata.urn,
                    start=sub_var["location"]["interval"]["start"]["value"],
                    end=sub_var["location"]["interval"]["end"]["value"],
                )
                if metadata.urn.startswith(
                    (
                        "urn:mavedb:00000047",
                        "urn:mavedb:00000048",
                        "urn:mavedb:00000053",
                        "urn:mavedb:00000058-a-1",
                    )
                ):
                    if tx_select_results is None:
                        raise ValueError  # these scoresets are all protein coding
                    seq = tx_select_results.sequence
                    sub_var["vrs_ref_allele_seq"] = seq[start_end[0] : start_end[1]]
                else:
                    seq = sub_var["location"]["sequence_id"]
                    sub_var["vrs_ref_allele_seq"] = sr.get_sequence(
                        seq, start_end[0], start_end[1]
                    )
        else:
            start_end = _format_start_end(
                metadata.urn,
                start=variant_list["location"]["interval"]["start"]["value"],
                end=variant_list["location"]["interval"]["end"]["value"],
            )
            if metadata.urn.startswith(
                (
                    "urn:mavedb:00000047",
                    "urn:mavedb:00000048",
                    "urn:mavedb:00000053",
                    "urn:mavedb:00000058-a-1",
                )
            ):
                if tx_select_results is None:
                    raise ValueError  # these scoresets are all protein coding
                seq = tx_select_results.sequence
                variant_list["vrs_ref_allele_seq"] = seq[start_end[0] : start_end[1]]
            else:
                seq = variant_list["location"]["sequence_id"]
                variant_list["vrs_ref_allele_seq"] = sr.get_sequence(
                    seq, start_end[0], start_end[1]
                )

        # Determine reference sequence
        if var.layer == AnnotationLayer.GENOMIC:
            if "members" in variant_list:
                acc = get_chromosome_identifier_from_vrs_id(
                    var.post_mapped_variants["members"][0]["location"]["sequence_id"]
                )
                if acc is None:
                    raise ValueError
                if acc.startswith("refseq:"):
                    acc = acc[7:]
            else:
                acc = get_chromosome_identifier_from_vrs_id(
                    var.post_mapped_variants["location"]["sequence_id"]
                )
                if acc is None:
                    raise ValueError
                if acc.startswith("refseq"):
                    acc = acc[7:]
        else:
            if tx_select_results is None:
                raise ValueError  # impossible by definition
            acc = tx_select_results.np

        # Add vrs_ref_allele_seq annotation and hgvs string to post-mapped variants
        variant_list = var.post_mapped_variants
        if "members" in variant_list:
            for sub_var in variant_list["members"]:
                sub_var["vrs_ref_allele_seq"] = sr.get_sequence(
                    sub_var["location"]["sequence_id"],
                    sub_var["location"]["interval"]["start"]["value"],
                    sub_var["location"]["interval"]["end"]["value"],
                )
                sub_var["hgvs"] = get_hgvs_string(sub_var, sr, acc)
        else:
            variant_list["vrs_ref_allele_seq"] = sr.get_sequence(
                variant_list["location"]["sequence_id"],
                variant_list["location"]["interval"]["start"]["value"],
                variant_list["location"]["interval"]["end"]["value"],
            )
            variant_list["hgvs"] = get_hgvs_string(variant_list, sr, acc)

    return vrs_results
