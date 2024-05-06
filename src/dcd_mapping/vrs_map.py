"""Map transcripts to VRS objects."""
import logging
from typing import List, Optional

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
    ScoreRow,
    ScoresetMetadata,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
    VrsMapping,
    VrsObject1_x,
)

__all__ = ["vrs_map"]


_logger = logging.getLogger(__name__)


class VrsMapError(Exception):
    """Raise in case of VRS mapping errors."""


def _create_hgvs_strings(
    alignment: AlignmentResult,
    raw_description: str,
    layer: AnnotationLayer,
    tx: Optional[TxSelectResult] = None,
) -> List[str]:
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
) -> Optional[VrsObject1_x]:
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
        return VrsMapping(
            mavedb_id=row.accession,
            pre_mapped_protein=[vrs_variation],
            post_mapped_protein=[vrs_variation],
            score=score,
        ).output_vrs_variations(AnnotationLayer.PROTEIN)
    layer = AnnotationLayer.PROTEIN
    hgvs_strings = _create_hgvs_strings(align_result, row.hgvs_pro, layer, transcript)
    return VrsMapping(
        mavedb_id=row.accession,
        score=score,
        pre_mapped_protein=_get_variation(
            hgvs_strings,
            layer,
            sequence_id,
            align_result,
            True,
        ),
        post_mapped_protein=_get_variation(
            hgvs_strings,
            layer,
            sequence_id,
            align_result,
            False,
            transcript.start,
        ),
    ).output_vrs_variations(AnnotationLayer.PROTEIN)


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
    records: List[ScoreRow],
    transcript: TxSelectResult,
    align_result: AlignmentResult,
) -> List[VrsObject1_x]:
    """Perform mapping on protein coding experiment results

    :param metadata: The metadata for a score set
    :param records: The list of MAVE variants in a given score set
    :param transcript: The transcript data for a score set
    :param align_results: The alignment data for a score set
    :return: A list of mappings
    """
    variations: List[VrsObject1_x] = []
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
        score = row.score
        hgvs_pro_mappings = _map_protein_coding_pro(
            row, score, align_result, psequence_id, transcript
        )
        if hgvs_pro_mappings:
            variations.append(hgvs_pro_mappings)
        if (
            row.hgvs_nt == "NA"
            or row.hgvs_nt in {"_wt", "_sy", "="}
            or len(row.hgvs_nt) == 3
        ):
            continue
        layer = AnnotationLayer.GENOMIC
        hgvs_strings = _create_hgvs_strings(align_result, row.hgvs_nt, layer)
        variations.append(
            VrsMapping(
                mavedb_id=row.accession,
                score=score,
                pre_mapped_genomic=_get_variation(
                    hgvs_strings,
                    layer,
                    gsequence_id,
                    align_result,
                    True,
                ),
                post_mapped_genomic=_get_variation(
                    hgvs_strings,
                    layer,
                    gsequence_id,
                    align_result,
                    False,
                ),
            ).output_vrs_variations(AnnotationLayer.GENOMIC)
        )
    return variations


def _map_regulatory_noncoding(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    align_result: AlignmentResult,
) -> List[VrsObject1_x]:
    """Perform mapping on noncoding/regulatory experiment results

    :param metadata: metadata for URN
    :param records: list of MAVE experiment result rows
    :param align_result: An AlignmentResult object for a score set
    :return: A list of VRS mappings
    """
    variations: List[VrsObject1_x] = []
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
        score = row.score
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
        variations.append(
            VrsMapping(
                pre_mapped_genomic=pre_map_allele,
                post_mapped_genomic=post_map_allele,
                mavedb_id=row.accession,
                score=score,
            ).output_vrs_variations(AnnotationLayer.GENOMIC)
        )
    return variations


def _get_variation(
    hgvs_strings: List[str],
    layer: AnnotationLayer,
    sequence_id: str,
    alignment: AlignmentResult,
    pre_map: bool,
    offset: int = 0,
) -> Optional[List[Allele]]:
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
    alleles: List[Allele] = []
    for hgvs_string in hgvs_strings:
        if (
            hgvs_string.endswith((".=", ")", "X")) or "?" in hgvs_string
        ):  # Invalid variant
            continue

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
                        alignment.query_subranges, alignment.hit_subranges
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
    records: List[ScoreRow],
    transcript: Optional[TxSelectResult] = None,
    silent: bool = True,
) -> Optional[List[VrsObject1_x]]:
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
