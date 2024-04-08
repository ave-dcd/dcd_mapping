"""Map transcripts to VRS objects.

Outstanding tasks/questions:
---------------------------
* Make sure typed digests vs full IDs vs plain digests are all being used correctly
* What is ``vrs_ref_allele_seq`` in output?
* Add basic transcript description where available
* Add HGVS expressions to alleles where available
"""
import functools
import logging
from typing import Dict, List, Optional, Union

import click
from Bio.Seq import Seq
from cool_seq_tool.schemas import AnnotationLayer, Strand
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs._internal.models import (
    Allele,
    Haplotype,
)
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
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
    VrsMappingResult,
)

__all__ = ["vrs_map"]


_logger = logging.getLogger(__name__)


class VrsMapError(Exception):
    """Raise in case of VRS mapping errors."""

def _create_hgvs_strings(
    alignment: AlignmentResult, raw_description: str, layer: AnnotationLayer
) -> List[str]:
    chrom_ac = get_chromosome_identifier(alignment.chrom)
    if '[' in raw_description:
        descr_list = list(set(raw_description[1:-1].split(";")))
        hgvs_strings = [f"{chrom_ac}:{layer.value}.{d}" for d in descr_list]
    else:
        descr_list = [raw_description]
        print(descr_list)
        hgvs_strings = [f"{chrom_ac}:{d}" for d in descr_list]
    return hgvs_strings

def _map_protein_coding_pro(
    row: ScoreRow,
    score: float,
    align_result: AlignmentResult,
    sequence: str,
    sequence_id: str,
    transcript: TxSelectResult,
) -> Optional[VrsMapping]:
    """Construct VRS object mapping for ``hgvs_pro`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row:
    :param score:
    :param align_result:
    :param sequence:
    :param sequence_id:
    :param transcript:
    :return: VRS mapping object if mapping succeeds
    """
    if row.hgvs_pro in {"_wt", "_sy", "p.="}:
        _logger.warning(
            f"Can't process Enrich2-style variant syntax {row.hgvs_nt} for {row.accession}"
        )
        return None
    if row.hgvs_pro.startswith("NP_009225.1:p."): # This is for experiment set 97
        vrs_variation = translate_hgvs_to_vrs(row.hgvs_pro)
        return VrsMapping(
            mavedb_id=row.accession,
            pre_mapped=vrs_variation,
            post_mapped=vrs_variation,
            score=score,
        )
    layer = AnnotationLayer.PROTEIN
    hgvs_strings = _create_hgvs_strings(align_result, row.hgvs_pro, layer)
    mapping = VrsMapping(
        mavedb_id=row.accession,
        score=score,
        pre_mapped=_get_variation(
            hgvs_strings,
            layer,
            sequence_id,
            sequence,
            align_result,
            True,
        ),
        post_mapped=_get_variation(
            hgvs_strings,
            layer,
            sequence_id,
            sequence,
            align_result,
            False,
            transcript.start,
        ),
    )
    return mapping

def _get_allele_sequence(allele: Allele) -> str:
    """Get sequence for Allele

    :param allele: VRS allele
    :return: sequence
    """
    sr = get_seqrepo()
    start = allele.location.start  # type: ignore
    end = allele.location.end  # type: ignore
    base = sr.sr[allele.location.sequenceReference.refgetAccession]  # type: ignore
    selection = base[start:end]
    return selection


def _map_protein_coding(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    transcript: TxSelectResult,
    align_result: AlignmentResult,
) -> VrsMappingResult:
    """Perform mapping on protein coding experiment results

    :param metadata:
    :param records:
    :param align_results:
    """
    variations = VrsMappingResult(variations=[])
    if metadata.target_sequence_type == TargetSequenceType.DNA:
        sequence = str(Seq(metadata.target_sequence).translate(table="1", stop_symbol=""))
    else:
        sequence = metadata.target_sequence

    sequence_id = f"SQ.{sha512t24u(sequence.encode('ascii'))}"
    alias_dict_list = [{'namespace': 'ga4gh', 'alias': sequence_id}]
    get_seqrepo().sr.store(sequence, nsaliases = alias_dict_list) # Add custom digest to SeqRepo

    for row in records:
        score = row.score

        # hgvs_pro
        hgvs_pro_mappings = _map_protein_coding_pro(
            row,
            score,
            align_result,
            sequence,
            sequence_id,
            transcript,
        )
        if hgvs_pro_mappings:
            variations.variations.append(hgvs_pro_mappings)
        if row.hgvs_nt == "NA":
            continue
        else:
            layer = AnnotationLayer.GENOMIC
            hgvs_strings = _create_hgvs_strings(align_result, row.hgvs_nt, layer)
            variations.variations.append(
                VrsMapping(
                    mavedb_id=row.accession,
                    score=score,
                    pre_mapped=_get_variation(
                        hgvs_strings,
                        layer,
                        sequence_id,
                        sequence,
                        align_result,
                        True,
                    ),
                    post_mapped=_get_variation(
                        hgvs_strings,
                        layer,
                        sequence_id,
                        sequence,
                        align_result,
                        False,
                    ),
                )
            )
    return variations

def _map_regulatory_noncoding(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    align_result: AlignmentResult,
) -> VrsMappingResult:
    """Perform mapping on noncoding/regulatory experiment results

    :param metadata: metadata for URN
    :param records: list of MAVE experiment result rows
    :param align_result:
    :return: TODO
    """
    variations = VrsMappingResult(variations=[])
    sequence_id = f"SQ.{sha512t24u(metadata.target_sequence.encode('ascii'))}"
    alias_dict_list = [{'namespace': 'ga4gh', 'alias': sequence_id}]
    get_seqrepo().sr.store(metadata.target_sequence, nsaliases = alias_dict_list) # Add custom digest to SeqRepo

    for row in records:
        if row.hgvs_nt in {"_wt", "_sy"}:
            _logger.warning(
                f"Can't process Enrich2-style variant syntax {row.hgvs_nt} for {metadata.urn}"
            )
            continue
        score = row.score
        hgvs_strings = _create_hgvs_strings(
            align_result, row.hgvs_nt[2:], AnnotationLayer.GENOMIC
        )
        pre_map_allele = _get_variation(
            hgvs_strings,
            AnnotationLayer.GENOMIC,
            sequence_id,
            metadata.target_sequence,
            align_result,
            pre_map=True,
            offset=0,
        )
        post_map_allele = _get_variation(
            hgvs_strings,
            AnnotationLayer.GENOMIC,
            sequence_id,
            metadata.target_sequence,
            align_result,
            pre_map=False,
            offset=0,
        )
        variations.variations.append(
            VrsMapping(
                pre_mapped=pre_map_allele,
                post_mapped=post_map_allele,
                mavedb_id=row.accession,
                score=score,
            )
        )
    return variations

def vrs_map(
    metadata: ScoresetMetadata,
    align_result: AlignmentResult,
    transcript: Optional[TxSelectResult],
    records: List[ScoreRow],
    silent: bool = True,
) -> Optional[VrsMappingResult]:
    """Given a description of a MAVE scoreset and an aligned transcript, generate
    the corresponding VRS objects.

    :param metadata: salient MAVE scoreset metadata
    :param transcript: output of transcript selection process
    :param records: scoreset records
    :param silent:
    :return: TODO
    """
    msg = f"Mapping {metadata.urn} to VRS..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)
    if metadata.target_gene_category == TargetType.PROTEIN_CODING and transcript:
        return _map_protein_coding(metadata, records, transcript, align_result)
    else:
        return _map_regulatory_noncoding(metadata, records, align_result)

def _get_variation(
    hgvs_strings: List[str],
    layer: AnnotationLayer,
    sequence_id: str,
    sequence: str,
    alignment: AlignmentResult,
    pre_map: bool,
    offset: int = 0,
) -> Union[Allele, List[Allele]]:
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

    :param hgvs_strings:
    :param layer: annotation layer
    :param sequence_id: target sequence digest eg ``"ga4gh:SQ.SQ.jUOcLPDjSqWFEo9kSOG8ITe1dr9QK3h6"``
    :param sequence: target sequence
    :param alignment:
    :param pre_map: if True, return object for pre mapping stage. Otherwise return for
        post-mapping.
    :param offset:
    :return:
    """
    if sequence_id.startswith("ga4gh:"):
        sequence_id = sequence_id[6:]
    alleles = []
    sequence_store = get_seqrepo()
    for hgvs_string in hgvs_strings:
        allele = translate_hgvs_to_vrs(hgvs_string)
        if "dup" in hgvs_string:
            allele.state.sequence = 2 * _get_allele_sequence(allele)  # type: ignore
        if pre_map:
            allele.location.sequenceReference.refgetAccession = sequence_id  # type: ignore
        else:
            if layer == AnnotationLayer.PROTEIN:
                allele.location.start += offset  # type: ignore
                allele.location.end += offset  # type: ignore
            else:
                start: int = allele.location.start  # type: ignore
                if (
                    len(alignment.query_subranges) == 1
                    and alignment.strand == Strand.POSITIVE
                ):
                    subrange_start = alignment.query_subranges[0].start
                    diff = start - subrange_start
                    diff2: int = allele.location.end - start  # type: ignore
                    allele.location.start = subrange_start + diff
                    allele.location.end = allele.location.start + diff2  # type: ignore
                else:
                    for query_subrange, hit_subrange in zip(
                        alignment.query_subranges, alignment.hit_subranges
                    ):
                        if start >= query_subrange.start and start < query_subrange.end:
                            query_subrange_start = query_subrange.start
                            hit_subrange_start = hit_subrange.start
                            break
                    else:
                        raise ValueError(
                            "Could not find hit subrange compatible with allele position"
                        )
                    diff = start - query_subrange_start
                    diff2: int = allele.location.end - start  # type: ignore
                    if alignment.strand == Strand.POSITIVE:  # positive orientation
                        allele.location.start = hit_subrange_start + diff
                        allele.location.end = allele.location.start + diff2  # type: ignore
                        if "dup" in hgvs_string:
                            allele.state.sequence = 2 * _get_allele_sequence(allele)  # type: ignore
                    else:
                        allele.location.start = hit_subrange_start - diff - diff2
                        allele.location.end = allele.location.start + diff2  # type: ignore
                        if "dup" in hgvs_string:
                            allele.state.sequence = 2 * _get_allele_sequence(allele)  # type: ignore
                        allele.state.sequence = str(
                            Seq(str(allele.state.sequence)).reverse_complement()
                        )
        if allele.state.sequence == "N" and layer != AnnotationLayer.PROTEIN:
            allele.state.sequence = _get_allele_sequence(allele)  # type: ignore
        if '=' in hgvs_string and layer == AnnotationLayer.PROTEIN:
            allele.state.sequence = _get_allele_sequence(allele)
        allele = normalize(allele, data_proxy=sequence_store)
        allele.id = ga4gh_identify(allele)
        alleles.append(allele)

    if len(alleles) == 1:
        return alleles[0]
    else:
        return alleles