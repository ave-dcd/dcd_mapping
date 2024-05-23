"""Map transcripts to VRS objects."""
import logging

import click
from Bio.Seq import Seq
from cool_seq_tool.schemas import AnnotationLayer, Strand
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs._internal.models import Allele, SequenceString
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
    VrsMapping,
)

__all__ = ["vrs_map"]


_logger = logging.getLogger(__name__)


class VrsMapError(Exception):
    """Raise in case of VRS mapping errors."""


def _create_hgvs_strings(
    alignment: AlignmentResult,
    raw_description: str,
    layer: AnnotationLayer,
    tx: TxSelectResult | None = None,
) -> list[str]:
    """Properly format MAVE variant strings
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
    return hgvs_strings


def _map_protein_coding_pro(
    row: ScoreRow,
    score: str,
    align_result: AlignmentResult,
    sequence_id: str,
    transcript: TxSelectResult,
) -> MappedScore | None:
    """Construct VRS object mapping for ``hgvs_pro`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row: A row of output from a MaveDB score set
    :param score: The score for a given row of output
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
            score=score,
            annotation_layer=AnnotationLayer.PROTEIN,
            pre_mapped=vrs_variation,
            post_mapped=vrs_variation,
            # TODO why were these defined as single-member arrays?
            # pre_mapped=[vrs_variation],
            # post_mapped=[vrs_variation]
        )
        return VrsMapping(
            mavedb_id=row.accession,
            pre_mapped_protein=[vrs_variation],
            post_mapped_protein=[vrs_variation],
            score=score,
        ).output_vrs_variations(AnnotationLayer.PROTEIN)
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
            score=score,
            annotation_layer=AnnotationLayer.PROTEIN,
            pre_mapped=pre_mapped_protein,
            post_mapped=post_mapped_protein,
        )
    return None
    # mapped = VrsMapping(
    #     mavedb_id=row.accession,
    #     score=score,
    #     pre_mapped_protein=_get_variation(
    #         hgvs_strings,
    #         AnnotationLayer.PROTEIN,
    #         sequence_id,
    #         align_result,
    #         True,
    #     ),
    #     post_mapped_protein=_get_variation(
    #         hgvs_strings,
    #         AnnotationLayer.PROTEIN,
    #         sequence_id,
    #         align_result,
    #         False,
    #         transcript.start,
    #     ),
    # )
    # if mapped.pre_mapped_protein and mapped.post_mapped_protein:
    #     return mapped.output_vrs_variations(AnnotationLayer.PROTEIN)
    # return None


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
    variations: list[MappedScore] = []
    if metadata.target_sequence_type == TargetSequenceType.DNA:
        sequence = str(
            Seq(metadata.target_sequence).translate(table="1", stop_symbol="")
        )
    else:
        sequence = metadata.target_sequence

    # Add custom digest to SeqRepo for both Protein and DNA Sequence
    psequence_id = store_sequence(sequence)
    gsequence_id = store_sequence(metadata.target_sequence)

    for row in records:
        hgvs_pro_mappings = _map_protein_coding_pro(
            row, row.score, align_result, psequence_id, transcript
        )
        if hgvs_pro_mappings:
            variations.append(hgvs_pro_mappings)
        if (
            row.hgvs_nt == "NA"
            or row.hgvs_nt in {"_wt", "_sy", "="}
            or len(row.hgvs_nt) == 3
        ):
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
        if pre_mapped_genomic is None or post_mapped_genomic is None:
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

        # variations.append(
        #     VrsMapping(
        #         mavedb_id=row.accession,
        #         score=score,
        #         pre_mapped_genomic=pre_mapped_genomic,
        #         post_mapped_genomic=post_mapped_genomic,
        #     ).output_vrs_variations(layer)
        # )
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
        if not pre_map_allele or not post_map_allele:
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
            # VrsMapping(
            #     pre_mapped_genomic=pre_map_allele,
            #     post_mapped_genomic=post_map_allele,
            #     mavedb_id=row.accession,
            #     score=score,
            # ).output_vrs_variations(AnnotationLayer.GENOMIC)
        )
    return variations


def _get_variation(
    hgvs_strings: list[str],
    layer: AnnotationLayer,
    sequence_id: str,
    alignment: AlignmentResult,
    pre_map: bool,
    offset: int = 0,
) -> list[Allele] | None:
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
    :return: A list of VRS Allele dictionaries for a list of MAVE variants
    """
    if sequence_id.startswith("ga4gh:"):
        sequence_id = sequence_id[6:]
    alleles: list[Allele] = []
    for hgvs_string in hgvs_strings:
        if hgvs_string.endswith((".=", ")", "X")):  # Invalid variant
            continue

        if "?" in hgvs_string:
            _logger.debug(
                "Substituting Xaa for ? in %s (sequence ID %s)",
                hgvs_string,
                sequence_id,
            )
            hgvs_string = hgvs_string.replace("?", "Xaa")

        # Generate VRS Allele structure. Set VA digests and SL digests to None
        allele = translate_hgvs_to_vrs(hgvs_string)
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
                                2 * _get_allele_sequence(allele)
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

        # Run ga4gh_identify to assign VA digest
        allele.id = ga4gh_identify(allele)
        alleles.append(allele)

    if not alleles:
        return None
    return alleles


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
    msg = f"Mapping {metadata.urn} to VRS..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)

    if metadata.urn == "urn:mavedb:00000072-a-1":
        msg = f"No RefSeq accession is available for {metadata.urn}."
        if not silent:
            click.echo(msg)
        _logger.warning(msg)
        return None

    if metadata.target_gene_category == TargetType.PROTEIN_CODING and transcript:
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
