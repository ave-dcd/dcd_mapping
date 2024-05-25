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
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import sha512t24u
from ga4gh.core._internal.models import Extension
from ga4gh.vrs._internal.models import (
    Allele,
    Expression,
    Haplotype,
    SequenceString,
    Syntax,
)

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
    MappedScore,
    ScoreAnnotation,
    ScoresetMetadata,
    TargetSequenceType,
    TxSelectResult,
    VrsMapping1_3,
    allele_to_vod,
    haplotype_to_haplotype_1_3,
)

_logger = logging.getLogger(__name__)


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


def _offset_allele_ref_seq(ss: str, start: int, end: int) -> tuple[int, int]:
    """Handle known edge cases in start and end coordinates for vrs_ref_allele_seq."""
    if ss.startswith("urn:mavedb:00000060-a-1"):
        _logger.warning(
            "urn:mavedb:00000060-a-1 reports the entire human reference sequence as the target sequence, but the start and end positions need to be manually offset by 289"
        )
        return (start + 289, end + 289)
    if ss.startswith("urn:mavedb:00000060-a-2"):
        _logger.warning(
            "urn:mavedb:00000060-a-2 reports the entire human reference sequence as the target sequence, but the start and end positions need to be manually offset by 331"
        )
        return (start + 331, end + 331)
    return (start, end)


def _get_vrs_ref_allele_seq(
    allele: Allele, metadata: ScoresetMetadata, tx_select_results: TxSelectResult | None
) -> Extension:
    start, end = _offset_allele_ref_seq(
        metadata.urn,
        allele.location.interval.start.value,  # type: ignore
        allele.location.interval.end.value,  # type: ignore
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
            # should be impossible - these scoresets are all protein coding
            raise ValueError
        seq = tx_select_results.sequence
        ref = seq[start:end]
    else:
        seq = f"ga4gh:{allele.location.sequenceReference.refgetAccession}"  # type: ignore
        sr = get_seqrepo()
        ref = sr.get_sequence(seq, start, end)
        if ref is None:
            raise ValueError
    return Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)


def _get_hgvs_string(allele: Allele, accession: str) -> tuple[str, Syntax]:
    """Return an HGVS string for a given VRS allele

    :param allele: A post-mapped VRS allele
    :param accession: A RefSeq accession
    :return An HGVS string and the Syntax value
    """
    if accession.startswith("NP"):
        syntax = Syntax.HGVS_P
        syntax_value = "p"
    else:
        syntax = Syntax.HGVS_G
        syntax_value = "g"
    start: int = allele.location.start  # type: ignore
    end: int = allele.location.end  # type: ignore

    dp = get_seqrepo()
    if start == end:
        ref = None
        aas = dp.get_sequence(accession, start - 1, start)
        aae = dp.get_sequence(accession, end, end + 1)
        end += 1
    else:
        ref = dp.get_sequence(accession, start, end)
        aas = dp.get_sequence(accession, start, start + 1)
        aae = dp.get_sequence(accession, end - 1, end)
        start += 1

    if syntax == Syntax.HGVS_P:
        ival = hgvs.location.Interval(
            start=hgvs.location.AAPosition(base=start, aa=aas),
            end=hgvs.location.AAPosition(base=end, aa=aae),
        )
    else:
        ival = hgvs.location.Interval(
            start=hgvs.location.SimplePosition(base=start),
            end=hgvs.location.SimplePosition(base=end),
        )
    alt: SequenceString = allele.state.sequence  # type: ignore

    edit = ""  # empty by default
    if alt == ref:
        edit = "="
    if ref and (2 * ref == alt or len(ref) == 1 and set(ref) == set(alt)):
        edit = "dup"
    if alt == "":
        edit = "del"

    if edit != "dup" or edit != "del" or edit != "=":
        edit = (
            hgvs.edit.AARefAlt(ref=ref, alt=alt)
            if syntax == Syntax.HGVS_P
            else hgvs.edit.NARefAlt(ref=ref, alt=alt)
        )

    if alt != ref:
        posedit = hgvs.posedit.PosEdit(pos=ival, edit=edit)
    else:
        posedit = (
            f"{seq3(ref)}{start!s}=" if syntax == Syntax.HGVS_P else f"{end!s}{ref}="
        )

    var = str(
        hgvs.sequencevariant.SequenceVariant(
            ac=accession, type=syntax_value, posedit=posedit
        )
    )
    if var.endswith("delins"):
        var = var.replace("delins", "del")
    return var, syntax


def _annotate_allele_mapping(
    mapped_score: MappedScore,
    tx_results: TxSelectResult | None,
    metadata: ScoresetMetadata,
) -> ScoreAnnotation:
    pre_mapped: Allele = mapped_score.pre_mapped
    post_mapped: Allele = mapped_score.post_mapped

    # get vrs_ref_allele_seq for pre-mapped variants
    pre_mapped.extensions = [_get_vrs_ref_allele_seq(post_mapped, metadata, tx_results)]

    # Determine reference sequence
    if mapped_score.annotation_layer == AnnotationLayer.GENOMIC:
        sequence_id = f"ga4gh:{mapped_score.post_mapped.location.sequenceReference.refgetAccession}"
        accession = get_chromosome_identifier_from_vrs_id(sequence_id)
        if accession is None:
            raise ValueError
        if accession.startswith("refseq:"):
            accession = accession[7:]
    else:
        if tx_results is None:
            raise ValueError  # impossible by definition
        accession = tx_results.np

    sr = get_seqrepo()
    loc = mapped_score.post_mapped.location
    sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"  # type: ignore
    ref = sr.get_sequence(sequence_id, loc.start, loc.end)  # TODO type issues???
    post_mapped.extensions = [
        Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)
    ]
    hgvs, syntax = _get_hgvs_string(post_mapped, accession)
    post_mapped.expressions = [Expression(syntax=syntax, value=hgvs)]

    pre_mapped_vod = allele_to_vod(pre_mapped)
    post_mapped_vod = allele_to_vod(post_mapped)

    return ScoreAnnotation(
        pre_mapped=pre_mapped_vod,
        post_mapped=post_mapped_vod,
        pre_mapped_2_0=pre_mapped,
        post_mapped_2_0=post_mapped,
        mavedb_id=mapped_score.accession_id,
        score=mapped_score.score,
    )


def _annotate_haplotype_mapping(
    mapping: MappedScore, tx_results: TxSelectResult | None, metadata: ScoresetMetadata
) -> ScoreAnnotation:
    pre_mapped: Haplotype = mapping.pre_mapped  # type: ignore
    post_mapped: Haplotype = mapping.post_mapped  # type: ignore
    allele: Allele
    # get vrs_ref_allele_seq for pre-mapped variants
    for allele in pre_mapped.members:
        allele.extensions = [_get_vrs_ref_allele_seq(allele, metadata, tx_results)]

    # Determine reference sequence
    if mapping.annotation_layer == AnnotationLayer.GENOMIC:
        sequence_id = (
            f"ga4gh:{post_mapped.members[0].location.sequenceReference.refgetAccession}"
        )
        accession = get_chromosome_identifier_from_vrs_id(sequence_id)
        if accession is None:
            raise ValueError
        if accession.startswith("refseq:"):
            accession = accession[7:]
    else:
        if tx_results is None:
            raise ValueError  # impossible by definition
        accession = tx_results.np

    sr = get_seqrepo()
    for allele in post_mapped.members:
        loc = allele.location
        sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"
        ref = sr.get_sequence(sequence_id, loc.start, loc.end)  # TODO type issues??
        allele.extensions = [
            Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)
        ]
        hgvs, syntax = _get_hgvs_string(allele, accession)
        allele.expressions = [Expression(syntax=syntax, value=hgvs)]

    pre_mapped_converted = haplotype_to_haplotype_1_3(pre_mapped)
    post_mapped_converted = haplotype_to_haplotype_1_3(post_mapped)

    return ScoreAnnotation(
        pre_mapped=pre_mapped_converted,
        post_mapped=post_mapped_converted,
        pre_mapped_2_0=pre_mapped,
        post_mapped_2_0=post_mapped,
        mavedb_id=mapping.accession_id,
        score=mapping.score,
    )


def annotate(
    mapped_scores: list[MappedScore],
    tx_results: TxSelectResult | None,
    metadata: ScoresetMetadata,
) -> list[ScoreAnnotation]:
    """Given a list of mappings, add additional contextual data:

    1. ``vrs_ref_allele_seq``: The sequence between the start and end positions
        indicated in the variant
    2. ``hgvs``: An HGVS string describing the variant (only included for post-mapped
        variants)
    3. ``transcript_accession``: A description of the MANE annotation of the transcript,
        if any  # < -- TODO I think we took this out

    ...and provide VRS 1.3-converted equivalents, too.

    :param vrs_results: in-progress variant mappings
    :param tx_select_results: transcript selection if available
    :param metadata: MaveDB scoreset metadata
    :return: annotated mappings objects
    """
    score_annotations = []
    for mapped_score in mapped_scores:
        if not mapped_score:
            msg = "#TODO: I would like to know when/why this happens"
            raise ValueError(msg)
        if isinstance(mapped_score.pre_mapped, Haplotype) and isinstance(
            mapped_score.post_mapped, Haplotype
        ):
            score_annotations.append(
                _annotate_haplotype_mapping(mapped_score, tx_results, metadata)
            )
        elif isinstance(mapped_score.pre_mapped, Allele) and isinstance(
            mapped_score.post_mapped, Allele
        ):
            score_annotations.append(
                _annotate_allele_mapping(mapped_score, tx_results, metadata)
            )
        else:
            ValueError("inconsistent variant structure")

    return score_annotations
