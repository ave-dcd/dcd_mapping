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
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import sha512t24u
from ga4gh.core._internal.models import Extension
from ga4gh.vrs._internal.models import (
    Allele,
    Expression,
    Haplotype,
    LiteralSequenceExpression,
    Syntax,
)

from dcd_mapping import vrs_v1_schemas
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


def _allele_to_v1_allele(allele: Allele) -> vrs_v1_schemas.Allele:
    """Convert VRS 2.0 allele to VRS 1.3 allele.

    :param allele: VRS 2.0a allele
    :return: equivalent VRS 1.3 allele
    """
    start = allele.location.start
    end = allele.location.end
    sequence_id = f"ga4gh:{allele.location.sequenceReference.refgetAccession}"  # type: ignore
    location_raw = f'{{"end":{{"type":"Number","value":{end}}},"sequence_id":"{sequence_id.split(".")[1]}","start":{{"type":"Number","value":{start}}},"type":"SequenceLocation"}}'
    location_id = sha512t24u(location_raw.encode("ascii"))
    sequence = "" if not allele.state.sequence else allele.state.sequence.root
    allele_raw = f'{{"location":"{location_id}","state":{{"sequence":"{sequence}","type":"LiteralSequenceExpression"}},"type":"Allele"}}'
    allele_id = sha512t24u(allele_raw.encode("ascii"))

    return vrs_v1_schemas.Allele(
        id=f"ga4gh:VA.{allele_id}",
        location=vrs_v1_schemas.SequenceLocation(
            id=location_id,
            sequence_id=sequence_id,
            interval=vrs_v1_schemas.SequenceInterval(
                start=vrs_v1_schemas.Number(value=start, type="Number"),
                end=vrs_v1_schemas.Number(value=end, type="Number"),
            ),
        ),
        state=vrs_v1_schemas.LiteralSequenceExpression(sequence=sequence),
    )


def _allele_to_vod(allele: Allele) -> vrs_v1_schemas.VariationDescriptor:
    """Convert VRS 2.0 allele to comparable VRSATILE VariationDescriptor.

    Some allele properties aren't available in the 1.3 allele, so we have to lift them
    up to the VariationDescriptor.
    """
    allele_v1 = _allele_to_v1_allele(allele)
    if allele.expressions:
        original_expression = allele.expressions[0]
        expressions = [
            vrs_v1_schemas.Expression(
                syntax=original_expression.syntax,
                value=original_expression.value,
                syntax_version=None,
            )
        ]
    else:
        expressions = []
    return vrs_v1_schemas.VariationDescriptor(
        id=allele_v1.id,
        variation=allele_v1,
        type="VariationDescriptor",
        expressions=expressions,
        vrs_ref_allele_seq=allele.extensions[0].value,
        extensions=[],
    )


def _haplotype_to_haplotype_1_3(haplotype: Haplotype) -> vrs_v1_schemas.Haplotype:
    """Convert VRS 2.0 Haplotype to VRS 1.3 Haplotype.

    :param haplotype: VRS 2.0 haplotype
    :return: VRS 1.3 haplotype that contains VRSATILE variation descriptors (not alleles)
    """
    members = []
    allele: Allele
    for allele in haplotype.members:  # type: ignore
        members.append(_allele_to_vod(allele))
    return vrs_v1_schemas.Haplotype(members=members)


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
) -> ScoreAnnotationWithLayer:
    """Perform annotations and create VRS 1.3 equivalents for allele mappings."""
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
    sequence_id = f"ga4gh:{loc.sequenceReference.refgetAccession}"
    ref = sr.get_sequence(sequence_id, loc.start, loc.end)
    post_mapped.extensions = [
        Extension(type="Extension", name="vrs_ref_allele_seq", value=ref)
    ]
    hgvs_string, syntax = _get_hgvs_string(post_mapped, accession)
    post_mapped.expressions = [Expression(syntax=syntax, value=hgvs_string)]

    pre_mapped_vod = _allele_to_vod(pre_mapped)
    post_mapped_vod = _allele_to_vod(post_mapped)

    return ScoreAnnotationWithLayer(
        pre_mapped=pre_mapped_vod,
        post_mapped=post_mapped_vod,
        pre_mapped_2_0=pre_mapped,
        post_mapped_2_0=post_mapped,
        mavedb_id=mapped_score.accession_id,
        score=float(mapped_score.score) if mapped_score.score else None,
        annotation_layer=mapped_score.annotation_layer,
    )


def _annotate_haplotype_mapping(
    mapping: MappedScore, tx_results: TxSelectResult | None, metadata: ScoresetMetadata
) -> ScoreAnnotationWithLayer:
    """Perform annotations and create VRS 1.3 equivalents for haplotype mappings."""
    pre_mapped: Haplotype = mapping.pre_mapped  # type: ignore
    post_mapped: Haplotype = mapping.post_mapped  # type: ignore
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

    pre_mapped_converted = _haplotype_to_haplotype_1_3(pre_mapped)
    post_mapped_converted = _haplotype_to_haplotype_1_3(post_mapped)

    return ScoreAnnotationWithLayer(
        pre_mapped=pre_mapped_converted,
        post_mapped=post_mapped_converted,
        pre_mapped_2_0=pre_mapped,
        post_mapped_2_0=post_mapped,
        mavedb_id=mapping.accession_id,
        score=float(mapping.score) if mapping.score is not None else None,
        annotation_layer=mapping.annotation_layer,
    )


def annotate(
    mapped_scores: list[MappedScore],
    tx_results: TxSelectResult | None,
    metadata: ScoresetMetadata,
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
    :return: annotated mappings objects
    """
    score_annotations = []
    for mapped_score in mapped_scores:
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
    include_vrs_2: bool = False,
    output_path: Path | None = None,
) -> Path:
    """Save mapping output for a score set in a JSON file

    :param urn: Score set accession
    :param mave_vrs_mappings: A dictionary of VrsObject1_x objects
    :param align_result: Alignment information for a score set
    :param tx_output: Transcript output for a score set
    :param include_vrs_2: if true, also include VRS 2.0 mappings
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
    mapped_scores: list[ScoreAnnotation] = []
    for m in mappings:
        if m.annotation_layer == preferred_layer:
            # drop annotation layer from mapping object
            mapped_scores.append(ScoreAnnotation(**m.model_dump()))

    output = ScoresetMapping(
        metadata=metadata,
        computed_reference_sequence=computed_reference_sequence,
        mapped_reference_sequence=mapped_reference_sequence,
        mapped_scores=mapped_scores,
    )

    if not include_vrs_2:
        for m in output.mapped_scores:
            m.pre_mapped_2_0 = None
            m.post_mapped_2_0 = None

    if not output_path:
        now = datetime.datetime.now(tz=datetime.timezone.utc).isoformat()
        output_path = LOCAL_STORE_PATH / f"{urn}_mapping_{now}.json"

    _logger.info("Saving mapping output to %s", output_path)
    with output_path.open("w") as file:
        json.dump(
            json.loads(output.model_dump_json(exclude_unset=True, exclude_none=True)),
            file,
            indent=4,
        )
    return output_path
