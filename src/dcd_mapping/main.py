"""Provide core MaveDB mapping methods."""
import logging
from pathlib import Path

import click

from dcd_mapping.align import AlignmentError, align
from dcd_mapping.annotate import annotate, save_mapped_output_json
from dcd_mapping.mavedb_data import get_scoreset_metadata, get_scoreset_records
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMetadata,
)
from dcd_mapping.transcripts import TxSelectError, select_transcript
from dcd_mapping.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


async def map_scoreset(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    silent: bool = True,
    output_path: Path | None = None,
) -> None:
    """Given information about a MAVE experiment, map to VRS and save output as JSON.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :param silent: if True, suppress console information output
    :param output_path: optional path to save output at
    """
    try:
        alignment_result = align(metadata, silent)
    except AlignmentError as e:
        _logger.error("Alignment failed for scoreset %s: %s", metadata.urn, e)
        raise e

    try:
        transcript = await select_transcript(
            metadata, records, alignment_result, silent
        )
    except TxSelectError as e:
        _logger.error("Transcript selection failed for scoreset %s", metadata.urn)
        raise e

    try:
        vrs_results = vrs_map(metadata, alignment_result, records, transcript, silent)
    except VrsMapError as e:
        _logger.error("VRS mapping failed for scoreset %s", metadata.urn)
        raise e
    if vrs_results is None:
        _logger.info("No mapping available for %s", metadata.urn)
        return

    vrs_results = annotate(transcript, vrs_results, metadata)
    save_mapped_output_json(
        metadata.urn, vrs_results, alignment_result, transcript, output_path
    )


async def map_scoreset_urn(
    urn: str, silent: bool = True, output_path: Path | None = None
) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param silent: if True, suppress console information output
    :param output_path: optional path to save output at
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, silent)
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        raise e
    await map_scoreset(metadata, records, silent, output_path)
