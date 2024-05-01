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
def cli(urn: str, debug: bool) -> None:
    """Get VRS mapping on preferred transcript for URN.

    For example:

    % dcd-map 'urn:mavedb:00000329-a-1'

    \f
    :param urn: scoreset URN
    :param debug: if True, enable debug logging
    """  # noqa: D301
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        filename="dcd-mapping.log",
        format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
        level=log_level,
        force=True,
    )
    _logger.debug("debug logging enabled")
    asyncio.run(map_scoreset_urn(urn, silent=False))


if __name__ == "__main__":
    cli()
