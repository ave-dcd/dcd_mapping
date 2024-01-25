"""Provide command-line interface for accessing mapping functions."""
import asyncio
import logging

import click

from dcd_mapping.main import map_scoreset_urn

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
    "--cache_align",
    "-c",
    is_flag=True,
    show_default=True,
    default=False,
    help="Enable caching for alignment results. Mostly useful for development/debugging.",
)
def cli(urn: str, debug: bool, cache_align: bool) -> None:
    """Get VRS mapping on preferred transcript for URN.

    For example:

    % dcd-map 'urn:mavedb:00000329-a-1'

    \f
    :param urn: scoreset URN
    :param debug: if True, enable debug logging
    :param cache_align: if True, save alignment output and reuse when available
    """  # noqa: D301
    if debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.basicConfig(
        filename="dcd-mapping.log",
        format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
        level=log_level,
        force=True,
    )
    _logger.debug("debug logging enabled")
    asyncio.run(map_scoreset_urn(urn, silent=False, cache_align=cache_align))


if __name__ == "__main__":
    cli()
