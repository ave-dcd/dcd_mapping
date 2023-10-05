"""Provide command-line interface for accessing mapping functions."""
import click

from .main import map_scoreset


@click.command(no_args_is_help=True)
@click.option("--urn", "-u", help="MaveDB scoreset URN")
def cli(urn: str) -> None:
    """Process user commands and call core `map_scoreset()` function.

    For example:

    % dcd-map --urn='urn:mavedb:00000329-a-1'

    \f
    :param urn: scoreset URN
    """  # noqa: D301
    map_scoreset(urn)


if __name__ == "__main__":
    cli()
