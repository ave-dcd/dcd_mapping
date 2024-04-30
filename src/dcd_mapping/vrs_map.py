"""Map transcripts to VRS objects.

Outstanding tasks/questions:
---------------------------
* Make sure typed digests vs full IDs vs plain digests are all being used correctly
* What is ``vrs_ref_allele_seq`` in output?
* Add basic transcript description where available
* Add HGVS expressions to alleles where available
"""
import logging
from typing import List, Optional, Union

import click
from Bio.Seq import Seq
from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.schemas import AnnotationLayer, Strand
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs._internal.models import Allele, SequenceString
from ga4gh.vrs.normalize import normalize

from dcd_mapping.lookup import (
    get_chromosome_identifier,
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
    VrsMappingResult,
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
    score: float,
    align_result: AlignmentResult,
    sequence: str,
    sequence_id: str,
    transcript: TxSelectResult,
    sr: SeqRepoAccess,
) -> Optional[VrsMapping]:
    """Construct VRS object mapping for ``hgvs_pro`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row: A row of output from a MaveDB score set
    :param score: The score for a given row of output
    :param align_result: The alignment data for a score set
    :param sequence: The target sequence for a score set
    :param sequence_id: The GA4GH accession for the provided sequence
    :param transcript: The transcript selection information for a score set
    :param sr: A SeqRepo object
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
    if row.hgvs_pro.startswith("NP_009225.1:p."):  # This is for experiment set 97
        vrs_variation = translate_hgvs_to_vrs(row.hgvs_pro)
        return VrsMapping(
            mavedb_id=row.accession,
            pre_mapped_protein=[vrs_variation],
            post_mapped_protein=[vrs_variation],
            score=score,
        ).output_vrs_variations(AnnotationLayer.PROTEIN, sr)
    layer = AnnotationLayer.PROTEIN
    hgvs_strings = _create_hgvs_strings(align_result, row.hgvs_pro, layer, transcript)
    return VrsMapping(
        mavedb_id=row.accession,
        score=score,
        pre_mapped_protein=_get_variation(
            hgvs_strings,
            layer,
            sequence_id,
            sequence,
            align_result,
            True,
        ),
        post_mapped_protein=_get_variation(
            hgvs_strings,
            layer,
            sequence_id,
            sequence,
            align_result,
            False,
            transcript.start,
        ),
    ).output_vrs_variations(AnnotationLayer.PROTEIN, sr)


def _get_allele_sequence(allele: Allele, sr: SeqRepoAccess) -> str:
    """Get sequence for Allele

    :param allele: VRS allele
    :param sr: SeqRepoAccess instance
    :return: sequence
    """
    start = allele.location.start
    end = allele.location.end
    base = sr.sr[f"ga4gh:{allele.location.sequenceReference.refgetAccession}"]
    return base[start:end]


def _map_protein_coding(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    transcript: TxSelectResult,
    align_result: AlignmentResult,
    sr: SeqRepoAccess,
) -> VrsMappingResult:
    """Perform mapping on protein coding experiment results

    :param metadata: The metadata for a score set
    :param records: The list of MAVE variants in a given score set
    :param transcript: The transcript data for a score set
    :param align_results: The alignment data for a score set
    :param sr: A SeqRepo object
    :return: A VrsMappingResult object
    """
    variations = VrsMappingResult(variations=[])
    if metadata.target_sequence_type == TargetSequenceType.DNA:
        sequence = str(
            Seq(metadata.target_sequence).translate(table="1", stop_symbol="")
        )
    else:
        sequence = metadata.target_sequence

    # Add custom digest to SeqRepo for both Protein and DNA Sequence
    psequence_id = f"SQ.{sha512t24u(sequence.encode('ascii'))}"
    alias_dict_list = [{"namespace": "ga4gh", "alias": psequence_id}]
    sr.sr.store(sequence, nsaliases=alias_dict_list)

    gsequence_id = f"SQ.{sha512t24u(metadata.target_sequence.encode('ascii'))}"
    alias_dict_list = [{"namespace": "ga4gh", "alias": gsequence_id}]
    sr.sr.store(metadata.target_sequence, nsaliases=alias_dict_list)

    for row in records:
        score = row.score
        # hgvs_pro
        hgvs_pro_mappings = _map_protein_coding_pro(
            row, score, align_result, sr, psequence_id, transcript, sr
        )
        if hgvs_pro_mappings:
            variations.variations.append(hgvs_pro_mappings)
        if (
            row.hgvs_nt == "NA"
            or row.hgvs_nt in {"_wt", "_sy", "="}
            or len(row.hgvs_nt) == 3
        ):
            continue
        layer = AnnotationLayer.GENOMIC
        hgvs_strings = _create_hgvs_strings(align_result, row.hgvs_nt, layer)
        variations.variations.append(
            VrsMapping(
                mavedb_id=row.accession,
                score=score,
                pre_mapped_genomic=_get_variation(
                    hgvs_strings,
                    layer,
                    gsequence_id,
                    sr,
                    align_result,
                    True,
                ),
                post_mapped_genomic=_get_variation(
                    hgvs_strings,
                    layer,
                    gsequence_id,
                    sr,
                    align_result,
                    False,
                ),
            ).output_vrs_variations(AnnotationLayer.GENOMIC, sr)
        )
    return variations


def _map_regulatory_noncoding(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    align_result: AlignmentResult,
    sr: SeqRepoAccess,
) -> VrsMappingResult:
    """Perform mapping on noncoding/regulatory experiment results

    :param metadata: metadata for URN
    :param records: list of MAVE experiment result rows
    :param align_result: An AlignmentResult object for a score set
    :param sr: A SeqRepo object
    :return: A VrsMappingResult object
    """
    variations = VrsMappingResult(variations=[])
    sequence_id = f"SQ.{sha512t24u(metadata.target_sequence.encode('ascii'))}"
    alias_dict_list = [{"namespace": "ga4gh", "alias": sequence_id}]
    sr.sr.store(
        metadata.target_sequence, nsaliases=alias_dict_list
    )  # Add custom digest to SeqRepo

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
            sr,
            align_result,
            pre_map=True,
            offset=0,
        )
        post_map_allele = _get_variation(
            hgvs_strings,
            AnnotationLayer.GENOMIC,
            sequence_id,
            sr,
            align_result,
            pre_map=False,
            offset=0,
        )
        variations.variations.append(
            VrsMapping(
                pre_mapped_genomic=pre_map_allele,
                post_mapped_genomic=post_map_allele,
                mavedb_id=row.accession,
                score=score,
            ).output_vrs_variations(AnnotationLayer.GENOMIC, sr)
        )
    return variations


def vrs_map(
    metadata: ScoresetMetadata,
    align_result: AlignmentResult,
    records: List[ScoreRow],
    sr: SeqRepoAccess,
    silent: bool = True,
    transcript: Optional[TxSelectResult] = None,
) -> Optional[VrsMappingResult]:
    """Given a description of a MAVE scoreset and an aligned transcript, generate
    the corresponding VRS objects.

    :param metadata: salient MAVE scoreset metadata
    :param align_result: output from the sequence alignment process
    :param records: scoreset records
    :param sr: A SeqRepo object
    :param silent: A boolean indicating whether output should be shown
    :param transcript: output of transcript selection process
    :return: A VrsMappingResult object
    """
    msg = f"Mapping {metadata.urn} to VRS..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)
    if metadata.target_gene_category == TargetType.PROTEIN_CODING and transcript:
        return _map_protein_coding(
            metadata=metadata,
            records=records,
            transcript=transcript,
            align_result=align_result,
            sr=sr,
        )
    return _map_regulatory_noncoding(
        metadata=metadata, records=records, align_result=align_result, sr=sr
    )


def _get_variation(
    hgvs_strings: List[str],
    layer: AnnotationLayer,
    sequence_id: str,
    sr: SeqRepoAccess,
    alignment: AlignmentResult,
    pre_map: bool,
    offset: int = 0,
) -> Union[List[dict]]:
    """Create variation (haplotype/allele).

    Outstanding questions:
    ---------------------
    * Make a class to shadow the seqrepo data proxy that handles custom sequences
    without adding them to the system seqrepo. As it currently stands, I wouldn't
    expect this code to complete successfully.
    * Do we really need to go through an HGVS string/the VRS translator? Can't we just
    build the object ourselves?
    * trim duplicate code
    * simply args

    :param hgvs_strings: The HGVS suffix that represents a variant
    :param layer: annotation layer
    :param sequence_id: target sequence digest eg ``"ga4gh:SQ.jUOcLPDjSqWFEo9kSOG8ITe1dr9QK3h6"``
    :param sequence: SeqRepo instance
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
    alleles = []
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
            allele.state.sequence = SequenceString(2 * _get_allele_sequence(allele, sr))
        if pre_map:
            allele.location.sequenceReference.refgetAccession = sequence_id
            if "dup" in hgvs_string:
                allele.state.sequence = SequenceString(
                    2 * _get_allele_sequence(allele, sr)
                )
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
                                2 * _get_allele_sequence(allele, sr)
                            )
                    else:
                        allele.location.start = hit_subrange.end - diff - diff2
                        allele.location.end = allele.location.start + diff2
                        if "dup" in hgvs_string:
                            allele.state.sequence = SequenceString(
                                2 * _get_allele_sequence(allele, sr)
                            )
                        temp_str = str(
                            Seq(str(allele.state.sequence.root)).reverse_complement()
                        )
                        allele.state.sequence = SequenceString(temp_str)
        if allele.state.sequence.root == "N" and layer == AnnotationLayer.GENOMIC:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele, sr))
        if "=" in hgvs_string and layer == AnnotationLayer.PROTEIN:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele, sr))
        allele = normalize(allele, data_proxy=sr)

        # Run ga4gh_identify to assign VA digest
        allele.id = ga4gh_identify(allele)
        alleles.append(allele)

    if not alleles:
        return None
    return alleles
