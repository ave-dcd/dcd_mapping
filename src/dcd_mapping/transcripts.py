"""Select best reference sequence.

Outstanding questions/confusion
-------------------------------
* ``urn:mavedb:00000097-n-1``: unable to find any matching transcripts
* Lots of scoresets (2023-) giving index errors in offset calculation
* Remove MANE sorting once upstream sorting is confirmed
"""
import logging
import re
from typing import List, Optional

from Bio.Data.CodonTable import IUPACData
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from cool_seq_tool.schemas import TranscriptPriority
from gene.database.database import click

from dcd_mapping.lookup import (
    get_chromosome_identifier,
    get_gene_symbol,
    get_mane_transcripts,
    get_protein_accession,
    get_seqrepo,
    get_sequence,
    get_transcripts,
    get_uniprot_sequence,
)
from dcd_mapping.schemas import (
    AlignmentResult,
    ManeDescription,
    ScoreRow,
    ScoresetMetadata,
    TargetSequenceType,
    TargetType,
    TranscriptDescription,
    TxSelectResult,
)

__all__ = ["select_transcript"]

_logger = logging.getLogger(__name__)


class TxSelectError(Exception):
    """Raise for transcript selection failure."""


async def _get_compatible_transcripts(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> List[List[str]]:
    """Acquire matching transcripts

    :param metadata: metadata for scoreset
    :param align_result: output of ``align()`` method
    :return: List of list of compatible transcripts
    """
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

    :param matching_transcripts: list of list of transcript accession IDs
    :return: list of transcripts shared by all sublists
    """
    common_transcripts_set = set(matching_transcripts[0])
    for sublist in matching_transcripts[1:]:
        common_transcripts_set.intersection_update(sublist)
    common_transcripts = list(common_transcripts_set)
    return common_transcripts


def _choose_best_mane_transcript(
    mane_transcripts: List[ManeDescription],
) -> Optional[ManeDescription]:
    """Choose best transcript (Select > Plus Clinical) given MANE status. This was
    originally a little longer but I think all we have to worry about is grabbing based
    on MANE status.

    TODO: this shouldn't be necessary anymore, we already sort them

    :param mane_transcripts: list of MANE transcript descriptions
    :return: best transcript
    """
    if not mane_transcripts:
        return None
    for transcript in mane_transcripts:
        if transcript.transcript_priority == TranscriptPriority.MANE_SELECT:
            return transcript
    for transcript in mane_transcripts:
        if transcript.transcript_priority == TranscriptPriority.MANE_PLUS_CLINICAL:
            return transcript
    return None


async def _get_longest_compatible_transcript(
    transcripts: List[str],
) -> Optional[TranscriptDescription]:
    """Get longest transcript from a list of compatible transcripts.

    I think there's a chance of some discord between UTA and Seqrepo and we might get
    KeyErrors here. If so, we should do a filter further up to drop any transcript
    accession IDs not recognized by SeqRepo.

    :param transcripts:
    :return:
    """
    transcripts.sort(key=lambda tx: len(get_seqrepo().sr[tx]))
    nm = transcripts[-1]
    np = await get_protein_accession(nm)
    if not np:
        return None
    return TranscriptDescription(
        refseq_nuc=nm,
        refseq_prot=np,
        transcript_priority=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )


def _get_protein_sequence(target_sequence: str) -> str:
    """Get protein sequence if necessary.

    It'd be nice if there was a more elegant way to check if the sequence was already a
    protein sequence.

    :param target_sequence: sequence set as baseline in MAVE experiment (might already
        be set to protein)
    :return: resulting protein sequence
    """
    if len(set(target_sequence)) > 4:
        protein_sequence = target_sequence
    else:
        protein_sequence = str(Seq(target_sequence).translate(table="1")).replace(
            "*", ""
        )
    return protein_sequence


async def _select_protein_reference(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> TxSelectResult:
    """Select preferred transcript for protein reference sequence

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: Best transcript and associated metadata
    :raise TxSelectError: if no matching MANE transcripts and unable to get UniProt ID/
    reference sequence
    """
    matching_transcripts = await _get_compatible_transcripts(metadata, align_result)
    if not matching_transcripts:
        common_transcripts = None
    else:
        common_transcripts = _reduce_compatible_transcripts(matching_transcripts)
    if not common_transcripts:
        if not metadata.target_uniprot_ref:
            raise TxSelectError(
                f"Unable to find matching transcripts for {metadata.urn}"
            )
        protein_sequence = get_uniprot_sequence(metadata.target_uniprot_ref.id)
        np_accession = metadata.target_uniprot_ref.id
        ref_sequence = get_uniprot_sequence(metadata.target_uniprot_ref.id)
        if not ref_sequence:
            raise ValueError(
                f"Unable to grab reference sequence from uniprot.org for {metadata.urn}"
            )
        nm_accession = None
        tx_mode = None
    else:
        mane_transcripts = get_mane_transcripts(common_transcripts)
        best_tx = _choose_best_mane_transcript(mane_transcripts)
        if not best_tx:
            best_tx = await _get_longest_compatible_transcript(common_transcripts)
        if not best_tx:
            raise TxSelectError
        ref_sequence = get_sequence(best_tx.refseq_prot)
        nm_accession = best_tx.refseq_nuc
        np_accession = best_tx.refseq_prot
        tx_mode = best_tx.transcript_priority

    protein_sequence = _get_protein_sequence(metadata.target_sequence)
    # TODO -- look at these two lines
    is_full_match = ref_sequence.find(protein_sequence) != -1
    start = ref_sequence.find(protein_sequence[:10])

    protein_mapping_info = TxSelectResult(
        nm=nm_accession,
        np=np_accession,
        start=start,
        is_full_match=is_full_match,
        sequence=protein_sequence,
        transcript_mode=tx_mode,
    )
    return protein_mapping_info


def _offset_target_sequence(metadata: ScoresetMetadata, records: List[ScoreRow]) -> int:
    """Find start location in target sequence (it's not always 0)

    :param metadata:
    :param records:
    :return: starting index position (may be 0)
    """
    if not isinstance(records[0].hgvs_pro, str) or records[0].hgvs_pro.startswith("NP"):
        return 0  # TODO explain?
    protein_change_list = [rec.hgvs_pro.lstrip("p.") for rec in records]

    # build table of parseable amino acids by reference location on target sequence
    amino_acids_by_position = {}
    for protein_change in protein_change_list:
        if protein_change == "_sy" or protein_change == "_wt":
            continue
        if ";" in protein_change:
            protein_changes = protein_change[1:-1].split(";")
        else:
            protein_changes = [protein_change]
        for change in protein_changes:
            aa = change[:3]
            if aa == "=" or change[-3:] not in IUPACData.protein_letters_3to1:
                continue
            if "=" in change:
                loc = change[3:-1]
            else:
                loc = change[3:-3]
            if loc not in amino_acids_by_position:
                loc = re.sub("[^0-9]", "", loc)
                if loc:
                    amino_acids_by_position[loc] = seq1(aa)

    err_locs = []
    protein_sequence = Seq(metadata.target_sequence).translate(table="1")
    for i in range(len(protein_sequence)):
        if (
            str(i) in amino_acids_by_position
            and amino_acids_by_position[str(i)] != protein_sequence[i - 1]
        ):
            err_locs.append(i)
    if len(err_locs) == 0:
        return 0

    amino_acids_by_position = {int(k): v for k, v in amino_acids_by_position.items()}
    amino_acids_by_position = sorted(amino_acids_by_position.items())
    amino_acids_by_position = dict(amino_acids_by_position)
    p0, p1, p2, p3, p4 = list(amino_acids_by_position.keys())[0:5]

    seq = ""
    for value in amino_acids_by_position.values():
        seq += value

    protein_sequence = _get_protein_sequence(metadata.target_sequence)
    offset = 0
    for i, base in enumerate(protein_sequence):
        if all(
            [
                base == amino_acids_by_position[p0],
                protein_sequence[i + p1 - p0] == amino_acids_by_position[p1],
                protein_sequence[i + p2 - p0] == amino_acids_by_position[p2],
                protein_sequence[i + p3 - p0] == amino_acids_by_position[p3],
                protein_sequence[i + p4 - p0] == amino_acids_by_position[p4],
            ]
        ):
            if i + 1 == min(amino_acids_by_position.keys()) or i + 2 == min(
                amino_acids_by_position.keys()
            ):
                offset = 0
            else:
                offset = i
    return offset


async def select_transcript(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    align_result: AlignmentResult,
    silent: bool = False,
) -> Optional[TxSelectResult]:
    """Select appropriate human reference sequence for scoreset.

    * Unnecessary for regulatory/other noncoding scoresets which report genomic
    variations.
    * For protein scoresets, identify a matching RefSeq protein reference sequence.

    :param metadata: Scoreset metadata from MaveDB
    :param records:
    :param align_result: alignment results
    :return: Transcript description (accession ID, offset, selected sequence, etc)
    """
    msg = f"Selecting reference sequence for {metadata.urn}..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)

    if metadata.urn in {"urn:mavedb:00000053-a-1", "urn:mavedb:00000053-a-2"}:
        # target sequence for these scoresets is missing codon
        return None

    if metadata.target_gene_category == TargetType.PROTEIN_CODING:
        transcript_reference = await _select_protein_reference(metadata, align_result)
        if (
            transcript_reference
            and metadata.target_sequence_type == TargetSequenceType.DNA
        ):
            offset = _offset_target_sequence(metadata, records)
            if offset:
                transcript_reference.start = offset
    else:
        # can't provide transcripts for regulatory/noncoding scoresets
        return None

    msg = "Reference selection complete."
    if not silent:
        click.echo(msg)
    return transcript_reference
