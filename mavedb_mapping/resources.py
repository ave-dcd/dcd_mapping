"""Manage external data resources.

This module is responsible for handling requests for MaveDB data, such as scoresets
or scoreset metadata. It should also instantiate any external resources needed for
tasks like transcript selection.

This isn't a priority, but eventually, methods that send data requests to
remote APIs should first check the local mavedb_mapping cache, and the
:py:module:`cache` module should be built out to support cache invalidation, remote
syncing, etc.
"""
import csv
import logging
from importlib import resources as impresources
from pathlib import Path
from typing import List, Set

import requests
from pydantic import ValidationError
from tqdm import tqdm

from mavedb_mapping.cache import LOCAL_STORE_PATH
from mavedb_mapping.schemas import ScoreRow, ScoresetMetadata

_logger = logging.getLogger("mavedb_mapping")


class ResourceAcquisitionError(Exception):
    """Raise when resource acquisition fails."""


def _http_download(url: str, out_path: Path, show_progress: bool = False) -> Path:
    """Download a file via HTTP.

    :param url: location of file to retrieve
    :param out_path: location to save file to
    :param show_progress: show TQDM progress bar if true
    :return: Path if download successful
    :raise requests.HTTPError: if request is unsuccessful
    """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get("content-length", 0))
        with open(out_path, "wb") as h:
            if show_progress:
                with tqdm(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    desc=f"Downloading {out_path.name}",
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


def fetch_all_scoreset_urns() -> Set[str]:
    """Fetch all scoreset URNs. Since species is annotated at the scoreset target level,
    we can't yet filter on anything like `homo sapien`.

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


def fetch_all_human_scoreset_urns() -> List[str]:
    """Fetch all human scoreset URNs. Pretty slow, shouldn't be used frequently because
    it requires requesting every single scoreset.

    :return: list of human scoreset URNs
    """
    scoreset_urns = fetch_all_scoreset_urns()
    human_scoresets: List[str] = []
    for urn in scoreset_urns:
        r = requests.get(f"https://api.mavedb.org/api/v1/score-sets/{urn}")
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError:
            _logger.info(f"Unable to retrieve scoreset data for URN {urn}")
            continue
        data = r.json()

        ref_maps = data.get("targetGene", {}).get("referenceMaps", [])
        if ref_maps:
            for ref_map in ref_maps:
                if ref_map.get("genome", {}).get("organismName", "") == "Homo sapiens":
                    human_scoresets.append(urn)
    return human_scoresets


def get_scoreset_metadata(scoreset_urn: str) -> ScoresetMetadata:
    """Acquire metadata for scoreset.

    Currently makes an API call every time. In the future, could be cached.

    :param scoreset_urn: URN for scoreset
    :return: Object containing salient metadata
    :raise ResourceAcquisitionError: if unable to acquire metadata
    """
    url = f"https://api.mavedb.org/api/v1/score-sets/{scoreset_urn}"
    r = requests.get(url)
    try:
        r.raise_for_status()
    except requests.HTTPError:
        _logger.error(f"Received HTTPError from {url}")
        raise ResourceAcquisitionError(f"Metadata for scoreset {scoreset_urn}")
    metadata = r.json()
    try:
        structured_data = ScoresetMetadata(
            urn=metadata["urn"],
            target_gene_name=metadata["targetGene"]["name"],
            target_gene_category=metadata["targetGene"]["category"],
            target_sequence=metadata["targetGene"]["wtSequence"]["sequence"],
            target_sequence_type=metadata["targetGene"]["wtSequence"]["sequenceType"],
        )
    except (KeyError, ValidationError) as e:
        _logger.error(
            f"Unable to extract metadata from API response for scoreset {scoreset_urn}: {e}"
        )
        raise ResourceAcquisitionError(f"Metadata for scoreset {scoreset_urn}")
    return structured_data


def get_scoreset_records(scoreset_urn: str) -> List[ScoreRow]:
    """Get scoreset records.

    Only hit the MaveDB API if unavailable locally. In the future, we should use
    caching utilities and function args to allow invalidation or force fetching from
    remote.

    :param scoreset_urn: URN for scoreset
    :return: Array of individual ScoreRow objects, containing information like protein
        changes and scores
    :raise ResourceAcquisitionError: if unable to fetch file
    """
    scores_csv = LOCAL_STORE_PATH / f"{scoreset_urn.replace(':', ' ')}_scores.csv"
    if not scores_csv.exists():
        url = f"https://api.mavedb.org/api/v1/score-sets/{scoreset_urn}/scores"
        try:
            _http_download(url, scores_csv)
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
    url: str = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit",
) -> Path:
    """Acquire reference genome file. This file shouldn't change, so no need
    to worry about cache handling once it's fetched.

    :param url: URL to fetch reference file from. By default, points to the USCS-hosted
        hg38 file in the 2bit file format.
    :return: path to acquired file
    :raise ResourceAcquisitionError: if unable to acquire file.
    """
    genome_file = LOCAL_STORE_PATH / "hg38.2bit"
    if not genome_file.exists():
        try:
            _http_download(url, genome_file, True)
        except requests.HTTPError:
            _logger.error(f"HTTPError when fetching reference genome file from {url}")
            raise ResourceAcquisitionError(f"Reference genome file at {url}")
    return genome_file


def get_mapping_tmp_dir() -> Path:
    """Acquire app-specific "tmp" directory. It's not actually temporary because it's
    manually maintained, but we need a slightly more durable file location than what the
    system tmp directory can provide.

    :return: path to temporary file directory
    """
    tmp: Path = impresources.files("mavedb_mapping") / "tmp"  # type: ignore
    tmp.mkdir(exist_ok=True)
    return tmp
