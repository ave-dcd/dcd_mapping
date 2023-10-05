"""Select best reference sequence."""
import logging
from typing import Dict, List

from Bio.Seq import Seq

from mavemap.lookup import (
    get_chromosome_identifier,
    get_gene_symbol,
    get_mane_transcripts,
    get_reference_sequence,
    get_transcripts,
)
from mavemap.schemas import (
    AlignmentResult,
    ManeData,
    ScoresetMetadata,
    TargetSequenceType,
    TargetType,
    TranscriptStatus,
)

_logger = logging.getLogger(__name__)


class TxSelectError(Exception):
    """Raise for transcript selection failure."""


async def _get_compatible_transcripts(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> List[List[str]]:
    """Acquire matching transcripts"""
    if align_result.chrom.startswith("chr"):
        aligned_chrom = align_result.chrom[3:]
    else:
        aligned_chrom = align_result.chrom
    chromosome = get_chromosome_identifier(aligned_chrom)
    gene_symbol = get_gene_symbol(metadata)
    if not gene_symbol:
        raise TxSelectError
    transcript_matches = []
    for hit_range in align_result.hit_subranges:
        matches_list = await get_transcripts(
            gene_symbol, chromosome, hit_range.start, hit_range.end
        )
        if matches_list:
            transcript_matches.append(matches_list)
    return transcript_matches


def _reduce_compatible_transcripts(matching_transcripts: List[List[str]]) -> List[str]:
    """Reduce list of list of transcripts to a list containing only entries present
    in each sublist

    :param tx_subranges_list: list of list of transcript accession IDs
    :return: list of transcripts shared by all sublists
    """
    common_transcripts_set = set(matching_transcripts[0])
    for sublist in matching_transcripts[1:]:
        common_transcripts_set.intersection_update(sublist)
    common_transcripts = list(common_transcripts_set)
    return common_transcripts


def _choose_best_transcript(mane_transcripts: List[ManeData], urn: str) -> ManeData:
    """Choose best transcript (Select > Plus Clinical) given MANE status

    Todo:
    ----
     * Handle case where it's empty (ie implement longest compatible logic)

    :param mane_transcripts: list of MANE transcript descriptions
    :param urn: scoreset ID for error message
    :return: best transcript
    """
    if len(mane_transcripts) == 2:
        if mane_transcripts[0].mane_status == TranscriptStatus.SELECT:
            return mane_transcripts[0]
        else:
            return mane_transcripts[1]
    elif len(mane_transcripts) == 1:
        return mane_transcripts[0]
    else:
        _logger.error(
            f"Unexpected number of MANE transcripts: {len(mane_transcripts)}, urn: {urn}"
        )
        raise TxSelectError


async def _select_protein_reference(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> Dict:
    """TODO -- return type is WIP

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: TODO
    """
    matching_transcripts = await _get_compatible_transcripts(metadata, align_result)
    common_transcripts = _reduce_compatible_transcripts(matching_transcripts)
    # TODO if this fails, look up transcripts on uniprot?

    mane_transcripts = get_mane_transcripts(common_transcripts)
    best_tx = _choose_best_transcript(mane_transcripts, metadata.urn)

    np, nm, status = best_tx.refseq_prot, best_tx.refseq_nuc, best_tx.mane_status

    # TODO there's a thing here about taking the sequence as-is if it contains
    # more than four unique chars, that seems off
    # Check specific chars used instead?
    if len(set(metadata.target_sequence)) > 4:
        protein_sequence = metadata.target_sequence
    else:
        protein_sequence = str(
            Seq(metadata.target_sequence).translate(table="1")
        ).replace("*", "")

    ref_sequence = get_reference_sequence(np)
    is_full_match = ref_sequence.find(protein_sequence)
    start = ref_sequence.find(protein_sequence[:10])  # TODO seems potentially sus?
    protein_mapping_info = {
        "nm": nm,
        "np": np,
        "start": start,
        "is_full_match": is_full_match,
        "transcript_mode": status,  # TODO ?????
    }
    return protein_mapping_info


async def select_reference(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> Dict:
    """Select appropriate human reference sequence for scoreset.

    Fairly trivial for regulatory/other noncoding scoresets which report genomic
    variations.
    For protein scoresets, identify a matching RefSeq protein reference sequence.
    More description here TODO.
    Return type unclear

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: TODO
    """
    if metadata.target_gene_category != TargetType.PROTEIN_CODING:
        raise ValueError  # TODO
    if metadata.target_sequence_type == TargetSequenceType.PROTEIN:
        return await _select_protein_reference(metadata, align_result)
    elif metadata.target_sequence_type == TargetSequenceType.DNA:
        return {}
    else:
        raise ValueError  # TODO
