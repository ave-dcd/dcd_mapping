"""Manage local data resources.

This module is responsible for handling requests for MaveDB data, such as scoresets
or scoreset metadata. It should also instantiate any external resources needed for
tasks like alignment. Finally, it also contains some methods and variables for
intermediary data (e.g. BLAT query/output files).

Much of this can/should be replaced by the ``mavetools`` library. (and/or
``wags-tails``.)
"""
import csv
import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set
from urllib.parse import urlparse

import click
import requests
from pydantic import ValidationError
from tqdm import tqdm

from dcd_mapping.schemas import ReferenceGenome, ScoreRow, ScoresetMetadata, UniProtRef

__all__ = [
    "get_cached_blat_output",
    "get_scoreset_urns",
    "get_human_urns",
    "get_raw_scoreset_metadata",
    "get_scoreset_records",
    "get_scoreset_metadata",
    "get_human_urns",
    "get_ref_genome_file",
    "get_mapping_tmp_dir",
]

_logger = logging.getLogger(__name__)


LOCAL_STORE_PATH = Path(
    os.environ.get(
        "MAVEDB_STORAGE_DIR", Path.home() / ".local" / "share" / "dcd-mapping"
    )
)
if not LOCAL_STORE_PATH.exists():
    LOCAL_STORE_PATH.mkdir(exist_ok=True, parents=True)


def get_cached_blat_output(urn: str) -> Optional[Path]:
    """Return cached BLAT output if available. Mostly useful for development/testing.

    :param urn: identifier for scoreset
    :return: path to BLAT output if exists
    """
    out_file = LOCAL_STORE_PATH / f"{urn}_blat_output.psl"
    if out_file.exists():
        return out_file
    else:
        return None


class ResourceAcquisitionError(Exception):
    """Raise when resource acquisition fails."""


def _http_download(url: str, out_path: Path, silent: bool = True) -> Path:
    """Download a file via HTTP.

    :param url: location of file to retrieve
    :param out_path: location to save file to
    :param silent: show TQDM progress bar if true
    :return: Path if download successful
    :raise requests.HTTPError: if request is unsuccessful
    """
    click.echo(f"Downloading {out_path.name} to {out_path.parents[0].absolute()}")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get("content-length", 0))
        with open(out_path, "wb") as h:
            if not silent:
                with tqdm(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    desc=out_path.name,
                    ncols=80,
                ) as progress_bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            h.write(chunk)
                            progress_bar.update(len(chunk))
            else:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        h.write(chunk)
    return out_path


def get_scoreset_urns() -> Set[str]:
    """Fetch all scoreset URNs. Since species is annotated at the scoreset target level,
    we can't yet filter on anything like `homo sapien` -- meaning this is fairly slow.

    :return: set of URN strings
    """
    r = requests.get("https://api.mavedb.org/api/v1/experiments/")
    r.raise_for_status()
    scoreset_urn_lists = [
        experiment["scoreSetUrns"]
        for experiment in r.json()
        if experiment.get("scoreSetUrns")
    ]
    return set([urn for urns in scoreset_urn_lists for urn in urns])


def _metadata_response_is_human(json_response: Dict) -> bool:
    """Check that response from scoreset metadata API refers to a human genome target.

    :param json_response: response from scoreset metadata API
    :return: True if contains a target tagged as ``"Homo sapiens"``
    """
    for target_gene in json_response.get("targetGenes", []):
        organism = (
            target_gene.get("targetSequence", {})
            .get("reference", {})
            .get("organismName")
        )
        if organism == "Homo sapiens":
            return True
    return False


def get_human_urns() -> List[str]:
    """Fetch all human scoreset URNs. Pretty slow, shouldn't be used frequently because
    it requires requesting every single scoreset.

    :return: list of human scoreset URNs
    """
    scoreset_urns = get_scoreset_urns()
    human_scoresets: List[str] = []
    for urn in scoreset_urns:
        r = requests.get(f"https://api.mavedb.org/api/v1/score-sets/{urn}")
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError:
            _logger.info(f"Unable to retrieve scoreset data for URN {urn}")
            continue
        data = r.json()

        if _metadata_response_is_human(data):
            human_scoresets.append(urn)
    return human_scoresets


def _get_uniprot_ref(scoreset_json: Dict[str, Any]) -> Optional[UniProtRef]:
    """Extract UniProt reference from scoreset metadata if available.

    :param scoreset_json: parsed JSON from scoresets API
    :return: UniProt ID if available
    """
    ext_ids = scoreset_json["targetGenes"][0].get("externalIdentifiers")
    if not ext_ids:
        return None
    for ext_id in ext_ids:
        if ext_id.get("identifier", {}).get("dbName") == "UniProt":
            return UniProtRef(
                id=f"uniprot:{ext_id['identifier']['identifier']}",
                offset=ext_id["offset"],
            )


def get_raw_scoreset_metadata(scoreset_urn: str) -> Dict:
    """Get raw (original JSON) metadata for scoreset.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param scoreset_urn: URN for scoreset
    :return: Complete JSON response for object
    :raise ResourceAcquisitionError: if HTTP request fails
    """
    metadata_file = LOCAL_STORE_PATH / f"{scoreset_urn}_metadata.json"
    if not metadata_file.exists():
        url = f"https://api.mavedb.org/api/v1/score-sets/{scoreset_urn}"
        r = requests.get(url)
        try:
            r.raise_for_status()
        except requests.HTTPError:
            _logger.error(f"Received HTTPError from {url}")
            raise ResourceAcquisitionError(f"Metadata for scoreset {scoreset_urn}")
        metadata = r.json()
        with open(metadata_file, "w") as f:
            json.dump(metadata, f)
    else:
        with open(metadata_file, "r") as f:
            metadata = json.load(f)
    return metadata


def get_scoreset_metadata(scoreset_urn: str) -> ScoresetMetadata:
    """Acquire metadata for scoreset.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param scoreset_urn: URN for scoreset
    :return: Object containing salient metadata
    :raise ResourceAcquisitionError: if unable to acquire metadata
    """
    metadata = get_raw_scoreset_metadata(scoreset_urn)
    if not _metadata_response_is_human(metadata):
        raise ResourceAcquisitionError(
            f"Experiment for {scoreset_urn} contains no human targets"
        )
    if len(metadata["targetGenes"]) > 1:
        raise ResourceAcquisitionError(
            f"Multiple target genes for {scoreset_urn} -- look into this."
        )
    gene = metadata["targetGenes"][0]
    try:
        structured_data = ScoresetMetadata(
            urn=metadata["urn"],
            target_gene_name=gene["name"],
            target_gene_category=gene["category"],
            target_sequence=gene["targetSequence"]["sequence"],
            target_sequence_type=gene["targetSequence"]["sequenceType"],
            target_reference_genome=gene["targetSequence"]["reference"]["shortName"],
            target_uniprot_ref=_get_uniprot_ref(metadata),
        )
    except (KeyError, ValidationError) as e:
        _logger.error(
            f"Unable to extract metadata from API response for scoreset {scoreset_urn}: {e}"
        )
        raise ResourceAcquisitionError(f"Metadata for scoreset {scoreset_urn}")

    return structured_data


def get_scoreset_records(scoreset_urn: str, silent: bool = True) -> List[ScoreRow]:
    """Get scoreset records.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param scoreset_urn: URN for scoreset
    :return: Array of individual ScoreRow objects, containing information like protein
        changes and scores
    :raise ResourceAcquisitionError: if unable to fetch file
    """
    scores_csv = LOCAL_STORE_PATH / f"{scoreset_urn}_scores.csv"
    # TODO use smarter/more flexible caching methods
    if not scores_csv.exists():
        url = f"https://api.mavedb.org/api/v1/score-sets/{scoreset_urn}/scores"
        try:
            _http_download(url, scores_csv, silent)
        except requests.HTTPError:
            _logger.error(f"HTTPError when fetching scores CSV from {url}")
            raise ResourceAcquisitionError(f"Scores CSV for scoreset {scoreset_urn}")

    scores_data: List[ScoreRow] = []
    with open(scores_csv, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            scores_data.append(ScoreRow(**row))
    return scores_data


def get_ref_genome_file(
    build: ReferenceGenome = ReferenceGenome.HG38, silent: bool = True
) -> Path:
    """Acquire reference genome file in 2bit format from UCSC.

    :param build: genome build to acquire
    :param silent: if True, suppress console output
    :return: path to acquired file
    :raise ResourceAcquisitionError: if unable to acquire file.
    """
    url = f"https://hgdownload.cse.ucsc.edu/goldenpath/{build.value.lower()}/bigZips/{build.value.lower()}.2bit"
    parsed_url = urlparse(url)
    genome_file = LOCAL_STORE_PATH / os.path.basename(parsed_url.path)
    # this file shouldn't change, so no need to think about more advanced caching
    if not genome_file.exists():
        try:
            _http_download(url, genome_file, silent)
        except requests.HTTPError:
            _logger.error(f"HTTPError when fetching reference genome file from {url}")
            raise ResourceAcquisitionError(f"Reference genome file at {url}")
    return genome_file


def get_mapping_tmp_dir() -> Path:
    """Acquire app-specific "tmp" directory. It's not actually temporary because it's
    manually maintained, but we need a slightly more durable file location than what the
    system tmp directory can provide. Used for storing small, consistently-named files
    like the BLAT query and results files.

    :return: path to temporary file directory
    """
    tmp = LOCAL_STORE_PATH / "tmp"
    tmp.mkdir(exist_ok=True)
    return tmp
