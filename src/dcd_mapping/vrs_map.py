"""Map transcripts to VRS objects."""

import logging
from itertools import cycle

import click
from Bio.Seq import Seq
from cool_seq_tool.schemas import AnnotationLayer, Strand
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs.models import (
    Allele,
    CisPhasedBlock,
    LiteralSequenceExpression,
    ReferenceLengthExpression,
    SequenceLocation,
    SequenceString,
)
from ga4gh.vrs.normalize import normalize

from dcd_mapping.lookup import (
    get_chromosome_identifier,
    get_seqrepo,
    translate_hgvs_to_vrs,
)
from dcd_mapping.schemas import (
    AlignmentResult,
    MappedScore,
    ScoreRow,
    ScoresetMetadata,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
)

__all__ = ["vrs_map", "VrsMapError"]


_logger = logging.getLogger(__name__)


class VrsMapError(Exception):
    """Raise in case of VRS mapping errors."""


def _hgvs_variant_is_valid(hgvs_string: str) -> bool:
    return not hgvs_string.endswith((".=", ")", "X"))


def _process_any_aa_code(hgvs_pro_string: str) -> str:
    """Substitute "Xaa" for "?" in variation expression.

    Some expressions seem to use the single-character "?" wildcard in the context of
    three-letter amino acid codes. This is weird, and the proper replacement is "Xaa".

    Note that we currently do NOT make any alterations to nucleotide strings that use
    weird apparently-wildcard characters like "X" -- we just treat them as invalid (see
    _hgvs_variant_is_valid()).

    :param hgvs_string: MAVE HGVS expression
    :return: processed variation (equivalent to input if no wildcard code found)
    """
    if "?" in hgvs_pro_string:
        _logger.debug("Substituting Xaa for ? in %s", hgvs_pro_string)
        hgvs_pro_string = hgvs_pro_string.replace("?", "Xaa")
    return hgvs_pro_string


def _create_hgvs_strings(
    alignment: AlignmentResult,
    raw_description: str,
    layer: AnnotationLayer,
    tx: TxSelectResult | None = None,
) -> list[str]:
    """Properly format MAVE variant strings

    * Add accession
    * Split up plural/'haplotype' variant expressions
    * Drop empty/nonexistent/non-computable variant expressions
    * Convert "?" -> "Xaa" (should only be for amino acid variation expressions)

    :param align_results: Alignment results for a score set
    :param raw_description: The variant list as expressed in MaveDB
    :param layer: The Annotation Layer (protein or genomic)
    :param tx: The transcript selection information for a score set
    :return A list of processed variants
    """
    if layer == AnnotationLayer.PROTEIN:
        if not tx:
            msg = "Can't get protein layer if no transcript selection results given"
            raise VrsMapError(msg)
        acc = tx.np
    else:
        acc = get_chromosome_identifier(alignment.chrom)
    if "[" in raw_description:
        descr_list = list(set(raw_description[3:-1].split(";")))
        hgvs_strings = [f"{acc}:{layer.value}.{d}" for d in descr_list]
    else:
        descr_list = [raw_description]
        hgvs_strings = [f"{acc}:{d}" for d in descr_list]
    hgvs_strings = list(filter(_hgvs_variant_is_valid, hgvs_strings))
    if layer == AnnotationLayer.PROTEIN:
        hgvs_strings = [_process_any_aa_code(s) for s in hgvs_strings]
    return hgvs_strings


def _map_protein_coding_pro(
    row: ScoreRow,
    align_result: AlignmentResult,
    sequence_id: str,
    transcript: TxSelectResult,
) -> MappedScore | None:
    """Construct VRS object mapping for ``hgvs_pro`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row: A row of output from a MaveDB score set
    :param align_result: The alignment data for a score set
    :param sequence: The target sequence for a score set
    :param sequence_id: The GA4GH accession for the provided sequence
    :param transcript: The transcript selection information for a score set
    :return: VRS mapping object if mapping succeeds
    """
    if (
        row.hgvs_pro in {"_wt", "_sy", "NA"}
        or "fs" in row.hgvs_pro
        or len(row.hgvs_pro) == 3
    ):
        _logger.warning(
            "Can't process variant syntax %s for %s", row.hgvs_pro, row.accession
        )
        return None
    # Special case for experiment set urn:mavedb:0000097
    if row.hgvs_pro.startswith("NP_009225.1:p."):
        vrs_variation = translate_hgvs_to_vrs(row.hgvs_pro)
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            annotation_layer=AnnotationLayer.PROTEIN,
            pre_mapped=vrs_variation,
            post_mapped=vrs_variation,
        )
    hgvs_strings = _create_hgvs_strings(
        align_result, row.hgvs_pro, AnnotationLayer.PROTEIN, transcript
    )
    pre_mapped_protein = _get_variation(
        hgvs_strings,
        AnnotationLayer.PROTEIN,
        sequence_id,
        align_result,
        True,
    )
    post_mapped_protein = _get_variation(
        hgvs_strings,
        AnnotationLayer.PROTEIN,
        sequence_id,
        align_result,
        False,
        transcript.start,
    )
    if pre_mapped_protein and post_mapped_protein:
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            annotation_layer=AnnotationLayer.PROTEIN,
            pre_mapped=pre_mapped_protein,
            post_mapped=post_mapped_protein,
        )
    return None


def _get_allele_sequence(allele: Allele) -> str:
    """Get sequence for Allele

    :param allele: VRS allele
    :return: sequence
    :raise ValueError: if sequence is none
    """
    dp = get_seqrepo()
    start = allele.location.start
    end = allele.location.end
    sequence = dp.get_sequence(
        f"ga4gh:{allele.location.sequenceReference.refgetAccession}", start, end
    )
    if sequence is None:
        raise ValueError
    return sequence


def store_sequence(sequence: str) -> str:
    """Store sequence in SeqRepo.

    :param sequence: raw sequence (ie nucleotides or amino acids)
    :return: sequence ID (sans prefix, which is ``"ga4gh"``)
    """
    sequence_id = f"SQ.{sha512t24u(sequence.encode('ascii'))}"
    alias_dict_list = [{"namespace": "ga4gh", "alias": sequence_id}]
    sr = get_seqrepo()
    sr.sr.store(sequence, alias_dict_list)
    return sequence_id


def _hgvs_nt_is_valid(hgvs_nt: str) -> bool:
    """Check for invalid or unavailable nucleotide MAVE-HGVS variation

    :param hgvs_nt: MAVE_HGVS nucleotide expression
    :return: True if expression appears populated and valid
    """
    return (
        (hgvs_nt != "NA")
        and (hgvs_nt not in {"_wt", "_sy", "="})
        and (len(hgvs_nt) != 3)
    )


def _map_protein_coding(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    transcript: TxSelectResult,
    align_result: AlignmentResult,
) -> list[MappedScore]:
    """Perform mapping on protein coding experiment results

    :param metadata: The metadata for a score set
    :param records: The list of MAVE variants in a given score set
    :param transcript: The transcript data for a score set
    :param align_results: The alignment data for a score set
    :return: A list of mappings
    """
    gsequence_id = ""
    psequence_id = ""
    if metadata.target_sequence_type == TargetSequenceType.DNA:
        if transcript:
            sequence = str(
                Seq(metadata.target_sequence).translate(table="1", stop_symbol="")
            )
            psequence_id = store_sequence(sequence)
        gsequence_id = store_sequence(metadata.target_sequence)
    else:
        sequence = metadata.target_sequence
        psequence_id = gsequence_id = store_sequence(sequence)

    variations: list[MappedScore] = []

    for row in records:
        if transcript:
            hgvs_pro_mappings = _map_protein_coding_pro(
                row, align_result, psequence_id, transcript
            )
            if hgvs_pro_mappings:
                variations.append(hgvs_pro_mappings)
        if not _hgvs_nt_is_valid(row.hgvs_nt):
            continue
        hgvs_strings = _create_hgvs_strings(
            align_result, row.hgvs_nt, AnnotationLayer.GENOMIC
        )
        pre_mapped_genomic = _get_variation(
            hgvs_strings,
            AnnotationLayer.GENOMIC,
            gsequence_id,
            align_result,
            True,
        )
        post_mapped_genomic = _get_variation(
            hgvs_strings,
            AnnotationLayer.GENOMIC,
            gsequence_id,
            align_result,
            False,
        )
        if pre_mapped_genomic is None and post_mapped_genomic is None:
            _logger.warning(
                "Encountered apparently invalid genomic variants in %s: %s",
                row.accession,
                row.hgvs_nt,
            )
            continue
        variations.append(
            MappedScore(
                accession_id=row.accession,
                score=row.score,
                annotation_layer=AnnotationLayer.GENOMIC,
                pre_mapped=pre_mapped_genomic,
                post_mapped=post_mapped_genomic,
            )
        )
    return variations


def _map_regulatory_noncoding(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    align_result: AlignmentResult,
) -> list[MappedScore]:
    """Perform mapping on noncoding/regulatory experiment results

    :param metadata: metadata for URN
    :param records: list of MAVE experiment result rows
    :param align_result: An AlignmentResult object for a score set
    :return: A list of VRS mappings
    """
    variations: list[MappedScore] = []
    sequence_id = store_sequence(metadata.target_sequence)

    for row in records:
        if (
            row.hgvs_nt in {"_wt", "_sy", "="}
            or "fs" in row.hgvs_nt
            or len(row.hgvs_nt) == 3
        ):
            _logger.warning(
                "Can't process variant syntax %s for %s", row.hgvs_nt, metadata.urn
            )
            continue
        hgvs_strings = _create_hgvs_strings(
            align_result, row.hgvs_nt, AnnotationLayer.GENOMIC
        )
        pre_map_allele = _get_variation(
            hgvs_strings,
            AnnotationLayer.GENOMIC,
            sequence_id,
            align_result,
            True,
            offset=0,
        )
        post_map_allele = _get_variation(
            hgvs_strings,
            AnnotationLayer.GENOMIC,
            sequence_id,
            align_result,
            False,
            offset=0,
        )
        if not pre_map_allele and not post_map_allele:
            msg = "Genomic variations missing"
            raise VrsMapError(msg)
        variations.append(
            MappedScore(
                accession_id=row.accession,
                annotation_layer=AnnotationLayer.GENOMIC,
                pre_mapped=pre_map_allele,
                post_mapped=post_map_allele,
                score=row.score,
            )
        )
    return variations


def _rle_to_lse(
    rle: ReferenceLengthExpression, location: SequenceLocation
) -> LiteralSequenceExpression:
    """Coerce ReferenceLengthExpression to LiteralSequenceExpression.

    RLEs are helpful for long repeating sequences but a) unnecessary here and b)
    create incompatibilities with some data extraction further down so to simplify,
    we'll just turn them into equivalent LiteralSequenceExpressions.
    """
    sr = get_seqrepo()
    sequence_id = location.sequenceReference.refgetAccession
    start: int = location.start
    end = start + rle.repeatSubunitLength
    subsequence = sr.get_sequence(f"ga4gh:{sequence_id}", start, end)
    c = cycle(subsequence)
    derived_sequence = ""
    for _ in range(rle.length):
        derived_sequence += next(c)
    return LiteralSequenceExpression(sequence=derived_sequence)


def _get_variation(
    hgvs_strings: list[str],
    layer: AnnotationLayer,
    sequence_id: str,
    alignment: AlignmentResult,
    pre_map: bool,
    offset: int = 0,
) -> Allele | CisPhasedBlock | None:
    """Create variation (allele).

    :param hgvs_strings: The HGVS suffix that represents a variant
    :param layer: annotation layer
    :param sequence_id: target sequence digest eg ``"ga4gh:SQ.jUOcLPDjSqWFEo9kSOG8ITe1dr9QK3h6"``
    :param alignment: The AlignmentResult object for a score set
    :param pre_map: if True, return object for pre mapping stage. Otherwise return for
        post-mapping.
    :param offset: The offset to adjust the start and end positions in allele. This
        parameter is used if the annotation layer is protein. For genomic variants, the
        offset is computed with respect to the alignment block.
    :return: an allele or cis-phased block of alleles
    """
    if sequence_id.startswith("ga4gh:"):
        sequence_id = sequence_id[6:]
    alleles: list[Allele] = []
    for hgvs_string in hgvs_strings:
        # Generate VRS Allele structure. Set VA digests and SL digests to None
        allele = translate_hgvs_to_vrs(hgvs_string)
        if allele is None:
            break
        allele.id = None
        allele.digest = None
        allele.location.id = None
        allele.location.digest = None

        if "dup" in hgvs_string:
            allele.state.sequence = SequenceString(2 * _get_allele_sequence(allele))
        if pre_map:
            allele.location.sequenceReference.refgetAccession = sequence_id
            if "dup" in hgvs_string:
                allele.state.sequence = SequenceString(2 * _get_allele_sequence(allele))
        else:
            if layer == AnnotationLayer.PROTEIN:
                allele.location.start += offset
                allele.location.end += offset
            else:
                start: int = allele.location.start
                if (
                    len(alignment.query_subranges) == 1
                    and alignment.strand == Strand.POSITIVE
                ):
                    subrange_start = alignment.query_subranges[0].start
                    diff = start - subrange_start
                    diff2: int = allele.location.end - start
                    allele.location.start = alignment.hit_subranges[0].start + diff
                    allele.location.end = allele.location.start + diff2
                else:
                    for query_subrange, hit_subrange in zip(  # noqa: B007  # TODO remove hit_subrange?
                        alignment.query_subranges, alignment.hit_subranges, strict=False
                    ):
                        if start >= query_subrange.start and start < query_subrange.end:
                            break
                    diff = start - query_subrange.start
                    diff2: int = allele.location.end - start
                    if alignment.strand == Strand.POSITIVE:  # positive orientation
                        allele.location.start = hit_subrange.start + diff
                        allele.location.end = allele.location.start + diff2
                        if "dup" in hgvs_string:
                            allele.state.sequence = SequenceString(
                                2 * _get_allele_sequence(allele)
                            )
                    else:
                        allele.location.start = hit_subrange.end - diff - diff2
                        allele.location.end = allele.location.start + diff2
                        if "dup" in hgvs_string:
                            allele.state.sequence = SequenceString(
                                str(
                                    Seq(
                                        2 * _get_allele_sequence(allele)
                                    ).reverse_complement()
                                )
                            )
                        temp_str = str(
                            Seq(str(allele.state.sequence.root)).reverse_complement()
                        )
                        allele.state.sequence = SequenceString(temp_str)
        if allele.state.sequence.root == "N" and layer == AnnotationLayer.GENOMIC:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele))
        if "=" in hgvs_string and layer == AnnotationLayer.PROTEIN:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele))
        allele = normalize(allele, data_proxy=get_seqrepo())
        if isinstance(allele.state, ReferenceLengthExpression):
            _logger.debug(
                "Coercing state for %s into LSE: %s",
                hgvs_string,
                allele.state.model_dump_json(),
            )
            allele.state = _rle_to_lse(allele.state, allele.location)

        # Run ga4gh_identify to assign VA digest
        allele.id = ga4gh_identify(allele)
        alleles.append(allele)

    if not alleles:
        return None
    if len(alleles) == 1:
        return alleles[0]
    cpb = CisPhasedBlock(members=alleles)
    cpb.id = ga4gh_identify(cpb)
    return cpb


def vrs_map(
    metadata: ScoresetMetadata,
    align_result: AlignmentResult,
    records: list[ScoreRow],
    transcript: TxSelectResult | None = None,
    silent: bool = True,
) -> list[MappedScore] | None:
    """Given a description of a MAVE scoreset and an aligned transcript, generate
    the corresponding VRS objects.

    :param metadata: salient MAVE scoreset metadata
    :param align_result: output from the sequence alignment process
    :param records: scoreset records
    :param transcript: output of transcript selection process
    :param silent: If true, suppress console output
    :return: A list of mapping results
    """
    if metadata.urn == "urn:mavedb:00000072-a-1":
        msg = f"No RefSeq accession is available for {metadata.urn}."
        if not silent:
            click.echo(msg)
        _logger.warning(msg)
        return None

    if metadata.target_gene_category == TargetType.PROTEIN_CODING:
        return _map_protein_coding(
            metadata,
            records,
            transcript,
            align_result,
        )
    return _map_regulatory_noncoding(
        metadata,
        records,
        align_result,
    )
