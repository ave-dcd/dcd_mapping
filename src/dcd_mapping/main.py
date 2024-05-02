"""Provide core MaveDB mapping methods."""
import logging
from typing import List

import click
from cool_seq_tool.schemas import AnnotationLayer

from dcd_mapping.align import AlignmentError, align
from dcd_mapping.lookup import get_chromosome_identifier_from_vrs_id, get_seqrepo
from dcd_mapping.mavedb_data import get_scoreset_metadata, get_scoreset_records
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMetadata,
    TxSelectResult,
    VrsMapping,
)
from dcd_mapping.transcripts import TxSelectError, select_transcript
from dcd_mapping.utils import get_hgvs_string
from dcd_mapping.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


def _format_start_end(ss: str, start: int, end: int) -> List[int]:
    """Format start and end coordinates for vrs_ref_allele_seq for known edge
    cases
    :param ss: score set
    :param start: start coordinate
    :param end: end coordinate
    :return A list of start and end coordinates
    """
    if ss.startswith("urn:mavedb:00000060-a-1"):
        # This score set set reports the entire human reference sequence as the
        # target sequence, but the positions in the score set occur with an offset
        # of 289
        return [start + 289, end + 289]
    if ss.startswith("urn:mavedb:00000060-a-2"):
        # This score set set reports the entire human reference sequence as the
        # target sequence, but the positions in the score set occur with an offset
        # of 331
        return [start + 331, end + 331]
    return [start, end]


def annotate(
    tx_select_results: TxSelectResult,
    vrs_results: List[VrsMapping],
    metadata: ScoresetMetadata,
) -> None:
    """TODO"""
    sr = get_seqrepo()
    for var in vrs_results:
        if not var:
            continue
        variant_list = var.pre_mapped_variants
        if "members" in variant_list:
            for sub_var in variant_list["members"]:
                start_end = _format_start_end(
                    metadata.urn,
                    start=sub_var["location"]["interval"]["start"]["value"],
                    end=sub_var["location"]["interval"]["end"]["value"],
                )
                if metadata.urn.startswith(
                    (
                        "urn:mavedb:00000047",
                        "urn:mavedb:00000048",
                        "urn:mavedb:00000053",
                        "urn:mavedb:00000058-a-1",
                    )
                ):
                    seq = tx_select_results.sequence
                    sub_var["vrs_ref_allele_seq"] = seq[start_end[0] : start_end[1]]
                else:
                    seq = sub_var["location"]["sequence_id"]
                    sub_var["vrs_ref_allele_seq"] = sr.get_sequence(
                        seq, start_end[0], start_end[1]
                    )
        else:
            start_end = _format_start_end(
                metadata.urn,
                start=variant_list["location"]["interval"]["start"]["value"],
                end=variant_list["location"]["interval"]["end"]["value"],
            )
            if metadata.urn.startswith(
                (
                    "urn:mavedb:00000047",
                    "urn:mavedb:00000048",
                    "urn:mavedb:00000053",
                    "urn:mavedb:00000058-a-1",
                )
            ):
                seq = tx_select_results.sequence
                variant_list["vrs_ref_allele_seq"] = seq[start_end[0] : start_end[1]]
            else:
                seq = variant_list["location"]["sequence_id"]
                variant_list["vrs_ref_allele_seq"] = sr.get_sequence(
                    seq, start_end[0], start_end[1]
                )

        # Determine reference sequence
        if var.layer == AnnotationLayer.GENOMIC:
            if "members" in variant_list:
                acc = get_chromosome_identifier_from_vrs_id(
                    var.post_mapped_variants["members"][0]["location"]["sequence_id"]
                )
                if acc.startswith("refseq:"):
                    acc = acc[7:]
            else:
                acc = get_chromosome_identifier_from_vrs_id(
                    var.post_mapped_variants["location"]["sequence_id"]
                )
                if acc.startswith("refseq"):
                    acc = acc[7:]
        else:
            acc = tx_select_results.np

        # Add vrs_ref_allele_seq annotation and hgvs string to post-mapped variants
        variant_list = var.post_mapped_variants
        if "members" in variant_list:
            for sub_var in variant_list["members"]:
                sub_var["vrs_ref_allele_seq"] = sr.get_sequence(
                    sub_var["location"]["sequence_id"],
                    sub_var["location"]["interval"]["start"]["value"],
                    sub_var["location"]["interval"]["end"]["value"],
                )
                sub_var["hgvs"] = get_hgvs_string(sub_var, sr, acc)
        else:
            variant_list["vrs_ref_allele_seq"] = sr.get_sequence(
                variant_list["location"]["sequence_id"],
                variant_list["location"]["interval"]["start"]["value"],
                variant_list["location"]["interval"]["end"]["value"],
            )
            variant_list["hgvs"] = get_hgvs_string(variant_list, sr, acc)


async def map_scoreset(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    silent: bool = True,
) -> None:
    """Given information about a MAVE experiment, map to VRS and save output as JSON.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :param silent: if True, suppress console output
    """
    try:
        alignment_result = align(metadata, silent)
    except AlignmentError as e:
        _logger.error("Alignment failed for scoreset %s: %s", metadata.urn, e)
        return

    try:
        transcript = await select_transcript(
            metadata, records, alignment_result, silent
        )
    except TxSelectError:
        _logger.error("Transcript selection failed for scoreset %s", metadata.urn)
        return

    try:
        vrs_results = vrs_map(metadata, alignment_result, records, transcript, silent)
    except VrsMapError:
        _logger.error("VRS mapping failed for scoreset %s", metadata.urn)
        return

    if transcript and vrs_results:
        annotate(transcript, vrs_results, metadata)


async def map_scoreset_urn(urn: str, silent: bool = True) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param silent: if True, suppress console output
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, silent)
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        return
    await map_scoreset(metadata, records, silent)
