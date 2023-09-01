"""Acquire MaveDB experiment and scoreset metadata."""
import csv
import logging
from pathlib import Path
from typing import List, Optional, Set

import requests
from pydantic import ValidationError

from mavedb_mapping.cache import LOCAL_STORE_PATH
from mavedb_mapping.schemas import ScoreRow, ScoresetMetadata

_logger = logging.getLogger("mavedb_mapping")


def get_all_scoreset_urns() -> Set[str]:
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


def get_all_human_scoreset_urns() -> List[str]:
    """Fetch all human scoreset URNs. Pretty slow, shouldn't be used frequently because
    it requires requesting every single scoreset.

    :return: list of human scoreset URNs
    """
    scoreset_urns = get_all_scoreset_urns()
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


def get_scoreset_metadata(scoreset_urn: str) -> Optional[ScoresetMetadata]:
    """Acquire metadata for scoreset.

    Runs an HTTP request every time. Ideally, we should allow the option to cache this,
    for users who want to work from a stable trove of data (e.g. for reproducibility) or
    to lessen pressure on the MaveDB API.

    :param scoreset_urn: URN for scoreset
    :return: Object containing salient metadata
    """
    r = requests.get(f"https://api.mavedb.org/api/v1/score-sets/{scoreset_urn}")
    r.raise_for_status()
    data = r.json()
    try:
        structured_data = ScoresetMetadata(
            urn=data["urn"],
            target_gene_name=data["targetGene"]["name"],
            target_gene_category=data["targetGene"]["category"],
            target_sequence=data["targetGene"]["wtSequence"]["sequence"],
            target_sequence_type=data["targetGene"]["wtSequence"]["sequenceType"],
        )
    except (KeyError, ValidationError):
        _logger.warning(
            f"Unable to extract metadata from API response for scoreset {scoreset_urn}"
        )
        return None
    return structured_data


def _get_scores_csv(scoreset_urn: str) -> Optional[Path]:
    """Acquire scoreset CSV file.

    Currently, if a local version is available, we'll use it. In the future, we should
    think about options to force fetch from MaveDB, or check if the cached version is
    invalidated.

    :param scoreset_urn: URN for scoreset
    :return: Path to downloaded CSV if acquisition succeeds
    """
    filename = f"{scoreset_urn.replace(':', ' ')}_scores.csv"
    outfile_path = LOCAL_STORE_PATH / filename
    if not outfile_path.exists():
        url = f"https://api.mavedb.org/api/v1/score-sets/{scoreset_urn}/scores"
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(outfile_path, "wb") as h:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        h.write(chunk)
    return outfile_path


def get_scores_data(scoreset_urn: str) -> List[ScoreRow]:
    """Get scoreset records.

    :param scoreset_urn: URN for scoreset
    :return: Array of individual ScoreRow objects, containing information like protein
        changes and scores
    """
    scores_csv = _get_scores_csv(scoreset_urn)
    if not scores_csv:
        raise Exception(f"Failed to acquire scores CSV for scoreset {scoreset_urn}")

    scores_data: List[ScoreRow] = []

    with open(scores_csv, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            scores_data.append(ScoreRow(**row))

    return scores_data
