"""Align MaveDB target sequences to a human reference genome."""
import subprocess
from pathlib import Path
from typing import List, Optional

from Bio.SearchIO import read as read_blat
from Bio.SearchIO._model import QueryResult
from Bio.SearchIO._model.hsp import HSP

from mavedb_mapping.resources import get_mapping_tmp_dir, get_ref_genome_file
from mavedb_mapping.schemas import (
    ScoresetMetadata,
    TargetSequenceType,
)


class AlignmentError(Exception):
    """Raise when errors encountered during alignment."""


def _build_query_file(scoreset_metadata: ScoresetMetadata) -> Path:
    """Construct BLAT query file.

    :param scoreset_metadata: MaveDB scoreset metadata object
    :return: Path to constructed file
    """
    query_file = get_mapping_tmp_dir() / "blat_query.fa"
    with open(query_file, "w") as f:
        f.write(">" + "query" + "\n")
        f.write(scoreset_metadata.target_sequence + "\n")
        f.close()
    return query_file


def _run_blat(scoreset_metadata: ScoresetMetadata) -> QueryResult:
    """Run a BLAT query and returns the output.

    How to handle the BLAT binary is a big TODO. There are some newish Python packages
    that wrap it/provide bindings -- might be best to look there.

    :param scoreset_metadata: object containing scoreset attributes
    :param ref_genome_path: location of reference genome file
    :return: BLAT query result
    :raise AlignmentError: if BLAT subprocess returns error code
    """
    query_file = _build_query_file(scoreset_metadata)
    reference_genome_file = get_ref_genome_file()
    min_score = len(scoreset_metadata.target_sequence) // 2  # minimum match 50%
    out_file = get_mapping_tmp_dir() / "blat_out.psl"

    if scoreset_metadata.target_sequence_type == TargetSequenceType.PROTEIN:
        command = f"blat {reference_genome_file} -q=prot -t=dnax -minScore={min_score} {query_file.absolute()} {out_file.absolute()}"
    else:
        # missing `-t=dnax`?
        command = f"blat {reference_genome_file} -minScore={min_score} {query_file.absolute()} {out_file.absolute()}"
    process = subprocess.run(
        command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
    )
    if process.returncode != 0:
        raise AlignmentError(
            f"BLAT process returned error code {process.returncode}: {command}"
        )

    # seems like maybe this hits an error resolved by adding -q=dnax sometimes?
    # investigate, refer to older code if it comes up
    output = read_blat("blat_out.psl", "blat-psl")
    return output


def _get_hit_starts(output: QueryResult, hit_index: int) -> List[int]:
    """Get hit starts from HSP object.

    :param output: BLAT result object
    :return: The starts of hit sequence.
    """
    hit_starts: List[int] = []
    for n in range(len(output[hit_index])):
        hit_starts.append(output[hit_index][n].hit_start)
    return hit_starts


def _get_hsp(output: QueryResult) -> HSP:
    """Obtain high-scoring pairs (HSP) for query sequence.

    (I am pretty sure this should be refactored)

    :param output: BLAT result object
    """
    hit_dict = {}
    for c in range(len(output)):
        max_hit = output[c][0].score
        for e in range(len(output[c])):
            if output[c][e].score > max_hit:
                max_hit = output[c][e].score
        hit_dict[c] = max_hit

    # Using top scoring hit
    hit = max(hit_dict, key=hit_dict.get)
    hit_starts = _get_hit_starts(output, hit)
    hsp = output[hit][hit_starts.index(max(hit_starts))]
    return hsp


def align(scoreset_metadata: ScoresetMetadata) -> Optional[HSP]:
    """Align target sequence to a reference genome.

    :param scoreset_metadata: object containing scoreset metadata
    :return: high-scoring pairs (HSP) object
    """
    try:
        blat_output = _run_blat(scoreset_metadata)
    except AlignmentError:
        return None
    return _get_hsp(blat_output)
