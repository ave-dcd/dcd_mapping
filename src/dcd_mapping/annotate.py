"""Annotate MaveDB score set metadata with mapped scores."""

import datetime
import json
import logging
from pathlib import Path

import hgvs.edit
import hgvs.location
import hgvs.parser
import hgvs.posedit
import hgvs.sequencevariant
from Bio.SeqUtils import seq3
from canonicaljson import encode_canonical_json
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import sha512t24u
from ga4gh.core.entity_models import Extension, Syntax
from ga4gh.core.identifiers import PrevVrsVersion, ga4gh_identify
from ga4gh.vrs.models import (
    Allele,
    CisPhasedBlock,
    Expression,
    LiteralSequenceExpression,
    SequenceString,
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
    MappedReferenceSequence,
    MappedScore,
    ScoreAnnotation,
    ScoreAnnotationWithLayer,
    ScoresetMapping,
    ScoresetMetadata,
    TargetSequenceType,
    TxSelectResult,
)

_logger = logging.getLogger(__name__)


def _get_vrs_1_3_ext(allele: Allele) -> Extension:
    return Extension(
        name="vrs_v1.3_id", value=ga4gh_identify(allele, as_version=PrevVrsVersion.V1_3)
    )


def _get_va_digest(allele: Allele) -> Extension:
    """Return the VA digest for a pre-mapped allele
    :param allele: A pre-mapped variant
    :return A VRS extension reporting the pre-mapped digest
    """
    return Extension(name="pre_mapped_id", value=allele.id)


def _is_valid_allele(allele: Allele, align_result: AlignmentResult) -> bool:
    """Check if a pre-mapped allele occurs within the alignment coverage
    :param allele: A pre-mapped allele
    :param align_result: Alignment data
    :return True if position occurs in coverage, False if not
    """
    return (
        align_result.query_range.start
        <= allele.location.start
        <= align_result.query_range.end
    )


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
    """Create `vrs_ref_allele_seq` property."""
    start, end = _offset_allele_ref_seq(
        metadata.urn, allele.location.start, allele.location.end
    )
    if (
        metadata.urn.startswith(
            (
                "urn:mavedb:00000047",
                "urn:mavedb:00000048",
                "urn:mavedb:00000053",
                "urn:mavedb:00000058-a-1",
            )
        )
        and tx_select_results is not None
    ):
        seq = tx_select_results.sequence
        ref = seq[start:end]
    else:
        seq = f"ga4gh:{allele.location.sequenceReference.refgetAccession}"  # type: ignore
        sr = get_seqrepo()
        ref = sr.get_sequence(seq, start, end)
        if ref is None:
            raise ValueError
    return SequenceString(root=ref)


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
    start: int = allele.location.start
    end: int = allele.location.end

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
    if isinstance(allele.state, LiteralSequenceExpression):
        alt = allele.state.sequence.root
    else:
        msg = (
            f"Unable to handle string for non-LSE based allele in {allele.model_dump()}"
        )
        raise NotImplementedError(msg)

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
    align_result: AlignmentResult,
) -> ScoreAnnotationWithLayer:
    """Perform annotations for allele mappings."""
    pre_mapped: Allele = mapped_score.pre_mapped
    post_mapped: Allele = mapped_score.post_mapped

    # get vrs_ref_allele_seq for pre-mapped variants
    pre_mapped.extensions = [
        _get_vrs_1_3_ext(pre_mapped),
    ]
    pre_mapped.location.sequence = _get_vrs_ref_allele_seq(
        pre_mapped, metadata, tx_results
    )

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
    sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"
    post_mapped.location.sequence = SequenceString(
        root=sr.get_sequence(sequence_id, loc.start, loc.end)
    )
    post_mapped.extensions = [
        _get_vrs_1_3_ext(post_mapped),
        _get_va_digest(pre_mapped),
    ]
    hgvs_string, syntax = _get_hgvs_string(post_mapped, accession)
    post_mapped.expressions = [Expression(syntax=syntax, value=hgvs_string)]

    namespace = metadata.urn
    val = mapped_score.accession_id.split("#")[1]

    # Check if post-mapped allele is valid
    if mapped_score.annotation_layer == AnnotationLayer.GENOMIC:
        post_mapped = (
            post_mapped if _is_valid_allele(pre_mapped, align_result) else None
        )

    # Remove extra digest attributes
    pre_mapped.digest = None
    if post_mapped:
        post_mapped.digest = None

    return ScoreAnnotationWithLayer(
        pre_mapped=pre_mapped,
        post_mapped=post_mapped,
        mavedb_id=f"{namespace}#{val}",
        score=float(mapped_score.score) if mapped_score.score else None,
        annotation_layer=mapped_score.annotation_layer,
    )


def _get_vrs_1_3_haplotype_id(cpb: CisPhasedBlock) -> str:
    allele_ids = sorted(
        [
            ga4gh_identify(a, as_version=PrevVrsVersion.V1_3).replace("ga4gh:VA.", "")
            for a in cpb.members
        ]
    )
    haplotype = {"members": allele_ids, "type": "Haplotype"}
    canonical_json = encode_canonical_json(haplotype)

    return f"ga4gh:VH.{sha512t24u(canonical_json)}"


def _annotate_cpb_mapping(
    mapping: MappedScore,
    tx_results: TxSelectResult | None,
    metadata: ScoresetMetadata,
    align_result: AlignmentResult,
) -> ScoreAnnotationWithLayer:
    """Perform annotations and create VRS 1.3 equivalents for CisPhasedBlock mappings."""
    pre_mapped: CisPhasedBlock = mapping.pre_mapped  # type: ignore
    post_mapped: CisPhasedBlock = mapping.post_mapped  # type: ignore
    # get vrs_ref_allele_seq for pre-mapped variants
    for allele in pre_mapped.members:
        allele.extensions = [
            _get_vrs_1_3_ext(allele),
        ]
        allele.location.sequence = _get_vrs_ref_allele_seq(allele, metadata, tx_results)
        allele.digest = None
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
    valid_post_mapped_alleles = []
    for post_mapped_allele, pre_mapped_allele in zip(
        post_mapped.members, pre_mapped.members, strict=True
    ):
        loc = post_mapped_allele.location
        sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"
        post_mapped_allele.location.sequence = SequenceString(
            root=sr.get_sequence(sequence_id, loc.start, loc.end)
        )
        post_mapped_allele.extensions = [
            _get_vrs_1_3_ext(post_mapped_allele),
            _get_va_digest(pre_mapped_allele),
        ]
        hgvs, syntax = _get_hgvs_string(post_mapped_allele, accession)
        post_mapped_allele.expressions = [Expression(syntax=syntax, value=hgvs)]
        if mapping.annotation_layer == AnnotationLayer.PROTEIN or _is_valid_allele(
            pre_mapped_allele, align_result
        ):
            valid_post_mapped_alleles.append(post_mapped_allele)
        post_mapped_allele.digest = None
    post_mapped.members = valid_post_mapped_alleles

    pre_mapped.extensions = [
        Extension(name="vrs_v1.3_id", value=_get_vrs_1_3_haplotype_id(pre_mapped))
    ]
    if len(post_mapped.members) >= 2:
        post_mapped.extensions = [
            Extension(name="vrs_v1.3_id", value=_get_vrs_1_3_haplotype_id(post_mapped)),
        ]
    elif len(post_mapped.members) == 1:
        post_mapped = post_mapped.members[0]
    else:
        post_mapped = None

    namespace = metadata.urn
    val = mapping.accession_id.split("#")[1]

    return ScoreAnnotationWithLayer(
        pre_mapped=pre_mapped,
        post_mapped=post_mapped,
        mavedb_id=f"{namespace}#{val}",
        score=float(mapping.score) if mapping.score is not None else None,
        annotation_layer=mapping.annotation_layer,
    )


def annotate(
    mapped_scores: list[MappedScore],
    tx_results: TxSelectResult | None,
    metadata: ScoresetMetadata,
    align_result: AlignmentResult,
) -> list[ScoreAnnotationWithLayer]:
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
    :param align_result: Alignment data
    :return: annotated mappings objects
    """
    score_annotations = []
    for mapped_score in mapped_scores:
        if isinstance(mapped_score.pre_mapped, CisPhasedBlock) and isinstance(
            mapped_score.post_mapped, CisPhasedBlock
        ):
            score_annotations.append(
                _annotate_cpb_mapping(mapped_score, tx_results, metadata, align_result)
            )
        elif isinstance(mapped_score.pre_mapped, Allele) and isinstance(
            mapped_score.post_mapped, Allele
        ):
            score_annotations.append(
                _annotate_allele_mapping(
                    mapped_score, tx_results, metadata, align_result
                )
            )
        else:
            ValueError("inconsistent variant structure")

    return score_annotations


def _get_computed_reference_sequence(
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


def _get_mapped_reference_sequence(
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


def _set_scoreset_layer(
    urn: str, mappings: list[ScoreAnnotationWithLayer]
) -> AnnotationLayer:
    """Many individual score results provide both genomic and protein variant
    expressions. If genomic expressions are available, that's what we'd like to use.
    This function tells us how to filter the final annotation objects.
    """
    if urn.startswith("urn:mavedb:00000097"):
        _logger.debug(
            "Manually selecting protein annotation for scores from urn:mavedb:00000097"
        )
        return AnnotationLayer.PROTEIN
    for mapping in mappings:
        if mapping.annotation_layer == AnnotationLayer.GENOMIC:
            return AnnotationLayer.GENOMIC
    return AnnotationLayer.PROTEIN


def save_mapped_output_json(
    urn: str,
    mappings: list[ScoreAnnotationWithLayer],
    align_result: AlignmentResult,
    tx_output: TxSelectResult | None,
    output_path: Path | None = None,
) -> Path:
    """Save mapping output for a score set in a JSON file

    :param urn: Score set accession
    :param mave_vrs_mappings: A dictionary of VrsObject1_x objects
    :param align_result: Alignment information for a score set
    :param tx_output: Transcript output for a score set
    :param output_path: specific location to save output to. Default to
        <dcd_mapping_data_dir>/urn:mavedb:00000XXX-X-X_mapping_<ISO8601 datetime>.json
    :return: output location
    """
    preferred_layer = _set_scoreset_layer(urn, mappings)
    metadata = get_raw_scoreset_metadata(urn)
    computed_reference_sequence = _get_computed_reference_sequence(
        urn, preferred_layer, tx_output
    )
    mapped_reference_sequence = _get_mapped_reference_sequence(
        preferred_layer, tx_output, align_result
    )
    mapped_scores: list[ScoreAnnotation] = [
        # drop annotation layer from mapping object
        ScoreAnnotation(**m.model_dump())
        for m in mappings
        if m.annotation_layer == preferred_layer
    ]

    output = ScoresetMapping(
        metadata=metadata,
        computed_reference_sequence=computed_reference_sequence,
        mapped_reference_sequence=mapped_reference_sequence,
        mapped_scores=mapped_scores,
    )

    if not output_path:
        now = datetime.datetime.now(tz=datetime.timezone.utc).isoformat()
        output_path = LOCAL_STORE_PATH / f"{urn}_mapping_{now}.json"

    _logger.info("Saving mapping output to %s", output_path)
    with output_path.open("w") as f:
        # temporarily using BaseModel.model_dump() -- should use .model_dump_json()
        # once fix to serializer is made in VRS-Python
        json.dump(output.model_dump(exclude_none=True), f, indent=4)

    return output_path
