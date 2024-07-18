"""Provide command-line interface for accessing mapping functions."""
import asyncio
import logging
from pathlib import Path

import click

from dcd_mapping.align import AlignmentError
from dcd_mapping.main import map_scoreset_urn
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.transcripts import TxSelectError
from dcd_mapping.vrs_map import VrsMapError

_logger = logging.getLogger(__name__)


@click.command(no_args_is_help=True)
@click.argument("urn", nargs=1)
@click.option(
    "--debug",
    "-d",
    is_flag=True,
    show_default=True,
    default=False,
    help="Enable debug logging",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Desired location at which output file should be saved",
)
@click.option(
    "--include_vrs_1_3",
    "-v",
    is_flag=True,
    default=False,
    help="Include VRS 1.3 mappings",
)
def cli(
    urn: str,
    debug: bool,
    output: Path | None,
    include_vrs_1_3: bool,
) -> None:
    """Get VRS mapping on preferred transcript for URN.

    For example:

    % dcd-map 'urn:mavedb:00000329-a-1'
    """
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        filename="dcd-mapping.log",
        format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
        level=log_level,
        force=True,
    )
    _logger.debug("debug logging enabled")
    try:
        asyncio.run(map_scoreset_urn(urn, output, include_vrs_1_3, silent=False))
    except (
        LookupError,
        AlignmentError,
        TxSelectError,
        VrsMapError,
        ResourceAcquisitionError,
    ):
        click.get_current_context().exit(1)


if __name__ == "__main__":
    cli()
