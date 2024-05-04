"""Provide core MaveDB mapping methods."""
import logging
from typing import List

import click

from dcd_mapping.align import AlignmentError, align
from dcd_mapping.annotate import annotate
from dcd_mapping.mavedb_data import get_scoreset_metadata, get_scoreset_records
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMetadata,
)
from dcd_mapping.transcripts import TxSelectError, select_transcript
from dcd_mapping.utils import save_mapped_output_json
from dcd_mapping.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


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
    if vrs_results is None:
        _logger.info("No mapping available for %s", metadata.urn)
        return

    vrs_results = annotate(transcript, vrs_results, metadata)
    save_mapped_output_json(
        ss=metadata.urn,
        mave_vrs_mappings=vrs_results,
        align_result=alignment_result,
        tx_output=transcript,
    )


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
