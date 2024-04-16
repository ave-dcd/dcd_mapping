"""Align MaveDB target sequences to a human reference genome."""
import logging
import os
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Any, Generator, List, Optional

from Bio.SearchIO import HSP
from Bio.SearchIO import read as read_blat
from Bio.SearchIO._model import Hit, QueryResult
from cool_seq_tool.schemas import Strand
from gene.database.database import click

from dcd_mapping.lookup import get_chromosome_identifier, get_gene_location
from dcd_mapping.resources import (
    LOCAL_STORE_PATH,
    get_mapping_tmp_dir,
    get_ref_genome_file,
)
from dcd_mapping.schemas import (
    AlignmentResult,
    GeneLocation,
    ScoresetMetadata,
    SequenceRange,
    TargetSequenceType,
)

__all__ = ["align"]

_logger = logging.getLogger(__name__)


class AlignmentError(Exception):
    """Raise when errors encountered during alignment."""


def _write_query_file(file: Path, lines: List[str]) -> None:
    """Write lines to query file. This method is broken out to enable easy mocking while
    testing.

    :param file: path to query file
    :param lines: list of lines to write (should be header and then sequence)
    """
    with file.open("w") as f:
        for line in lines:
            f.write(f"{line}\n")


def _build_query_file(scoreset_metadata: ScoresetMetadata) -> Generator[Path, Any, Any]:
    """Construct BLAT query file.

    :param scoreset_metadata: MaveDB scoreset metadata object
    :return: Yielded Path to constructed file. Deletes file once complete.
    """
    query_file = (
        get_mapping_tmp_dir() / f"blat_query_{scoreset_metadata.urn}_{uuid.uuid1()}.fa"
    )
    _logger.debug("Writing BLAT query to %s", query_file)
    lines = [">query", scoreset_metadata.target_sequence]
    _write_query_file(query_file, lines)
    yield query_file
    query_file.unlink()


def _run_blat(
    target_args: str, query_file: Path, out_file: str, silent: bool
) -> subprocess.CompletedProcess:
    """Execute BLAT binary with relevant params.

    Currently, we rely on a system-installed BLAT binary accessible in the containing
    environment's PATH, or under env var ``BLAT_BIN_PATH``. This is sort of awkward and
    it'd be nice to make use of some direct bindings or better packaging if that's possible.

    * Perhaps `gget`? https://pachterlab.github.io/gget/en/blat.html
    * ``PxBlat``? https://github.com/ylab-hi/pxblat

    :param target_args: target params eg ``"-q=prot -t=dnax"`` (can be empty)
    :param query_file: path to query FASTA file
    :param out_file: path-like string to output fill (could be "/dev/stdout")
    :param silent: if True, suppress all console output
    :return: process result
    """
    reference_genome_file = get_ref_genome_file(silent=silent)
    bin_name = os.environ["BLAT_BIN_PATH"] if "BLAT_BIN_PATH" in os.environ else "blat"  # noqa: SIM401
    command = f"{bin_name} {reference_genome_file} {target_args} -minScore=20 {query_file} {out_file}"
    _logger.debug("Running BLAT command: %s", command)
    result = subprocess.run(  # noqa: UP022
        command,
        shell=True,  # noqa: S602
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _logger.debug("BLAT command finished with result %s", result.returncode)
    if result.returncode != 0:
        msg = f"BLAT process returned error code {result.returncode}: {target_args} {query_file} {out_file}"
        raise AlignmentError(msg)
    return result


def _get_cached_blat_output(metadata: ScoresetMetadata, silent: bool) -> QueryResult:
    """Get a BLAT output object for the given scoreset -- either reusing a previously-
    run query output file, or making a new one if unavailable.

    This method is broken out from ``_get_blat_output`` because it was getting pretty
    messy to handle differing file management logic between the two of them. However,
    query/command generation logic should be shared.

    :param metadata:
    :param silent:
    :return:
    """
    out_file = LOCAL_STORE_PATH / f"{metadata.urn}_blat_output.psl"
    if out_file.exists():
        return read_blat(out_file.absolute(), "blat-psl")
    query_file = next(_build_query_file(metadata))
    if metadata.target_sequence_type == TargetSequenceType.PROTEIN:
        target_args = "-q=prot -t=dnax"
    else:
        target_args = ""
    _run_blat(target_args, query_file, str(out_file.absolute()), silent)
    try:
        output = read_blat(out_file.absolute(), "blat-psl")
    except ValueError:
        target_args = "-q=dnax -t=dnax"
        _run_blat(target_args, query_file, str(out_file.absolute()), silent)
        try:
            output = read_blat(out_file.absolute(), "blat-psl")
        except ValueError as e:
            msg = f"Unable to get valid BLAT response for {metadata.urn}"
            raise AlignmentError(msg) from e
    return output


def _write_blat_output_tempfile(result: subprocess.CompletedProcess) -> str:
    """Create temp BLAT output file. Not immediately deleted, but should eventually
    be cleared by the OS.

    :param result: BLAT process result object
    :return: path-like string representing file location
    """
    raw_output = result.stdout.split(b"Loaded")[0]
    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.write(raw_output)
    return tmp.name


def _get_blat_output(metadata: ScoresetMetadata, silent: bool) -> QueryResult:
    """Run a BLAT query and returns a path to the output object.

    If unable to produce a valid query the first time, then try a query using ``dnax``
    bases.

    :param scoreset_metadata: object containing scoreset attributes
    :param silent: suppress BLAT command output
    :return: BLAT query result
    :raise AlignmentError: if BLAT subprocess returns error code
    """
    query_file = next(_build_query_file(metadata))
    if metadata.target_sequence_type == TargetSequenceType.PROTEIN:
        target_args = "-q=prot -t=dnax"
    else:
        target_args = ""
    process_result = _run_blat(target_args, query_file, "/dev/stdout", silent)
    out_file = _write_blat_output_tempfile(process_result)

    try:
        output = read_blat(out_file, "blat-psl")
    except ValueError:
        target_args = "-q=dnax -t=dnax"
        process_result = _run_blat(target_args, query_file, "/dev/stdout", silent)
        out_file = _write_blat_output_tempfile(process_result)
        try:
            output = read_blat(out_file, "blat-psl")
        except ValueError as e:
            msg = f"Unable to run successful BLAT on {metadata.urn}"
            raise AlignmentError(msg) from e

    return output


def _get_best_hit(output: QueryResult, urn: str, chromosome: Optional[str]) -> Hit:
    """Get best hit from BLAT output.

    First, try to return hit corresponding to expected chromosome taken from scoreset
    metadata. If chromosome doesn't match any of the outputs or is unavailable, take
    the hit with the single highest-scoring HSP.

    :param output: BLAT output
    :param urn: scoreset URN to use in error messages
    :param chromosome: refseq chromosome ID, e.g. ``"NC_000001.11"``
    :return: best Hit
    :raise AlignmentError: if unable to get hits from output
    """
    if chromosome:
        if chromosome.startswith("refseq"):
            chromosome = chromosome[7:]

        for hit in output:
            hit_chr = hit.id
            if hit_chr.startswith("chr"):
                hit_chr = hit_chr[3:]
            hit_chr_ac = get_chromosome_identifier(hit_chr)
            if hit_chr_ac == chromosome:
                return hit
        else:
            if list(output):
                hit_chrs = [h.id for h in output]
                _logger.warning(
                    "Failed to match hit chromosomes during alignment. URN: %s, expected chromosome: %s, hit chromosomes: %s",
                    urn,
                    chromosome,
                    hit_chrs,
                )

    best_score = 0
    best_score_hit = None
    for hit in output:
        best_local_score = max(hit, key=lambda i: i.score).score
        if best_local_score > best_score:
            best_score = best_local_score
            best_score_hit = hit

    if best_score_hit is None:
        _logger.error("Couldn't get hits from %s -- check BLAT output.", urn)
        raise AlignmentError

    return best_score_hit


def _get_best_hsp(hit: Hit, urn: str, gene_location: Optional[GeneLocation]) -> HSP:
    """Retrieve preferred HSP from BLAT Hit object.

    If gene location data is available, prefer the HSP with the least distance
    between the start of the hit and the start coordinate of the gene. Otherwise,
    take the HSP with the highest score value.

    :param hit: hit object from BLAT result
    :param urn: scoreset identifier for use in error messages
    :param gene_location: location data acquired by normalizing scoreset metadata
    :return: Preferred HSP object
    :raise AlignmentError: if hit object appears to be empty (should be impossible)
    """
    best_hsp = None
    if gene_location and gene_location.start is not None:
        best_hsp = min(hit, key=lambda hsp: abs(hsp.hit_start - gene_location.start))
    else:
        best_hsp = max(hit, key=lambda hsp: hsp.score)
    if best_hsp is None:
        _logger.error(
            "Unable to get best HSP from hit -- this should be impossible? urn: %s, hit: %s",
            urn,
            hit,
        )
        raise AlignmentError
    return best_hsp


def _get_best_match(output: QueryResult, metadata: ScoresetMetadata) -> AlignmentResult:
    """Obtain best high-scoring pairs (HSP) object for query sequence.

    :param metadata: scoreset metadata
    :param output: BLAT result object
    :return: alignment result ??
    """
    location = get_gene_location(metadata)
    chromosome = location.chromosome if location else None
    best_hit = _get_best_hit(output, metadata.urn, chromosome)
    best_hsp = _get_best_hsp(best_hit, metadata.urn, location)

    strand = Strand.POSITIVE if best_hsp[0].query_strand == 1 else Strand.NEGATIVE
    coverage = 100 * (best_hsp.query_end - best_hsp.query_start) / output.seq_len
    identity = best_hsp.ident_pct
    chrom = best_hsp.hit_id

    query_subranges = []
    hit_subranges = []
    for fragment in best_hsp:
        query_subranges.append(
            SequenceRange(start=fragment.query_start, end=fragment.query_end)
        )
        hit_subranges.append(
            SequenceRange(start=fragment.hit_start, end=fragment.hit_end)
        )

    return AlignmentResult(
        chrom=chrom,
        strand=strand,
        ident_pct=identity,
        coverage=coverage,
        query_range=SequenceRange(start=best_hsp.query_start, end=best_hsp.query_end),
        query_subranges=query_subranges,
        hit_range=SequenceRange(start=best_hsp.hit_start, end=best_hsp.hit_end),
        hit_subranges=hit_subranges,
    )


def align(
    scoreset_metadata: ScoresetMetadata, silent: bool = True, use_cached: bool = False
) -> AlignmentResult:
    """Align target sequence to a reference genome.

    :param scoreset_metadata: object containing scoreset metadata
    :param silent: suppress BLAT process output if true
    :param use_cached: make use of permanent mapping storage for intermediary files rather
        than rebuilding new output and storing in tmp directory. Mostly useful for
        development/testing.
    :return: data wrapper containing alignment results
    """
    msg = f"Performing alignment for {scoreset_metadata.urn}..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)

    if use_cached:
        blat_output = _get_cached_blat_output(scoreset_metadata, silent)
    else:
        blat_output = _get_blat_output(scoreset_metadata, silent)
    match = _get_best_match(blat_output, scoreset_metadata)

    msg = "Alignment complete."
    if not silent:
        click.echo(msg)
    _logger.info(msg)

    return match
