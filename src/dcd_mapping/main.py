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


def _emit_info(msg: str, silent: bool) -> None:
    if not silent:
        click.echo(msg)
    _logger.info(msg)


async def map_scoreset(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    output_path: Path | None = None,
    include_vrs_2: bool = False,
    silent: bool = True,
) -> None:
    """Given information about a MAVE experiment, map to VRS and save output as JSON.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :param output_path: optional path to save output at
    :param include_vrs_2: if true, include VRS 2.0 mappings in output JSON
    :param silent: if True, suppress console information output
    """
    _emit_info(f"Performing alignment for {metadata.urn}...", silent)
    try:
        alignment_result = align(metadata, silent)
    except AlignmentError as e:
        _logger.error("Alignment failed for scoreset %s: %s", metadata.urn, e)
        raise e
    _emit_info("Alignment complete.", silent)

    _emit_info("Selecting reference sequence...", silent)
    try:
        transcript = await select_transcript(metadata, records, alignment_result)
    except TxSelectError as e:
        _logger.error("Transcript selection failed for scoreset %s", metadata.urn)
        raise e
    _emit_info("Reference selection complete.", silent)

    _emit_info("Mapping to VRS...", silent)
    try:
        vrs_results = vrs_map(metadata, alignment_result, records, transcript, silent)
    except VrsMapError as e:
        _logger.error("VRS mapping failed for scoreset %s", metadata.urn)
        raise e
    if vrs_results is None:
        _logger.info("No mapping available for %s", metadata.urn)
        return
    _emit_info("VRS mapping complete.", silent)

    _emit_info("Annotating metadata and saving to file...", silent)
    vrs_results = annotate(vrs_results, transcript, metadata)
    final_output = save_mapped_output_json(
        metadata.urn,
        vrs_results,
        alignment_result,
        transcript,
        include_vrs_2,
        output_path,
    )
    _emit_info(f"Annotated scores saved to: {final_output}.", silent)


async def map_scoreset_urn(
    urn: str,
    output_path: Path | None = None,
    include_vrs_2: bool = False,
    silent: bool = True,
) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param output_path: optional path to save output at
    :param include_vrs_2: if true, include VRS 2.0 mappings in output JSON
    :param silent: if True, suppress console information output
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, silent)
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        raise e
    await map_scoreset(metadata, records, output_path, include_vrs_2, silent)
