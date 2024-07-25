"""Handle requests for MaveDB data, such as scoresets or scoreset metadata.

Much of this can/should be replaced by the ``mavetools`` library? (and/or ``wags-tails``.)
"""

import csv
import json
import logging
import tempfile
import zipfile
from pathlib import Path
from typing import Any

import requests
from pydantic import ValidationError

from dcd_mapping.resource_utils import (
    LOCAL_STORE_PATH,
    ResourceAcquisitionError,
    http_download,
)
from dcd_mapping.schemas import ScoreRow, ScoresetMetadata, UniProtRef

__all__ = [
    "get_scoreset_urns",
    "get_human_urns",
    "get_raw_scoreset_metadata",
    "get_scoreset_records",
    "get_scoreset_metadata",
    "get_human_urns",
]

_logger = logging.getLogger(__name__)


def get_scoreset_urns() -> set[str]:
    """Fetch all scoreset URNs. Since species is annotated at the scoreset target level,
    we can't yet filter on anything like `homo sapien` -- meaning this is fairly slow.

    :return: set of URN strings
    """
    r = requests.get("https://api.mavedb.org/api/v1/experiments/", timeout=30)
    r.raise_for_status()
    scoreset_urn_lists = [
        experiment["scoreSetUrns"]
        for experiment in r.json()
        if experiment.get("scoreSetUrns")
    ]
    return {urn for urns in scoreset_urn_lists for urn in urns}


def _metadata_response_is_human(json_response: dict) -> bool:
    """Check that response from scoreset metadata API refers to a human genome target.

    :param json_response: response from scoreset metadata API
    :return: True if contains a target tagged as ``"Homo sapiens"``
    """
    for target_gene in json_response.get("targetGenes", []):
        organism = (
            target_gene.get("targetSequence", {})
            .get("taxonomy", {})
            .get("organismName")
        )
        if organism == "Homo sapiens":
            return True
    return False


def get_human_urns() -> list[str]:
    """Fetch all human scoreset URNs. Pretty slow, shouldn't be used frequently because
    it requires requesting every single scoreset.

    :return: list of human scoreset URNs
    """
    scoreset_urns = get_scoreset_urns()
    human_scoresets: list[str] = []
    for urn in scoreset_urns:
        r = requests.get(f"https://api.mavedb.org/api/v1/score-sets/{urn}", timeout=30)
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError:
            _logger.info("Unable to retrieve scoreset data for URN %s", urn)
            continue
        data = r.json()

        if _metadata_response_is_human(data):
            human_scoresets.append(urn)
    return human_scoresets


def _get_uniprot_ref(scoreset_json: dict[str, Any]) -> UniProtRef | None:
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
    return None


def get_raw_scoreset_metadata(
    scoreset_urn: str, dcd_mapping_dir: Path | None = None
) -> dict:
    """Get raw (original JSON) metadata for scoreset.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param scoreset_urn: URN for scoreset
    :param dcd_mapping_dir: optionally declare location to save metadata to.
    :return: Complete JSON response for object
    :raise ResourceAcquisitionError: if HTTP request fails
    """
    if not dcd_mapping_dir:
        dcd_mapping_dir = LOCAL_STORE_PATH
    metadata_file = dcd_mapping_dir / f"{scoreset_urn}_metadata.json"
    if not metadata_file.exists():
        url = f"https://api.mavedb.org/api/v1/score-sets/{scoreset_urn}"
        r = requests.get(url, timeout=30)
        try:
            r.raise_for_status()
        except requests.HTTPError as e:
            msg = f"Received HTTPError from {url} for scoreset {scoreset_urn}"
            _logger.error(msg)
            raise ResourceAcquisitionError(msg) from e
        metadata = r.json()
        with metadata_file.open("w") as f:
            json.dump(metadata, f)
    else:
        with metadata_file.open() as f:
            metadata = json.load(f)
    return metadata


def get_scoreset_metadata(
    scoreset_urn: str, dcd_mapping_dir: Path | None = None
) -> ScoresetMetadata:
    """Acquire metadata for scoreset.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param scoreset_urn: URN for scoreset
    :param dcd_mapping_dir: optionally declare location to save metadata to.
    :return: Object containing salient metadata
    :raise ResourceAcquisitionError: if unable to acquire metadata
    """
    metadata = get_raw_scoreset_metadata(scoreset_urn, dcd_mapping_dir)

    if not _metadata_response_is_human(metadata):
        msg = f"Experiment for {scoreset_urn} contains no human targets"
        raise ResourceAcquisitionError(msg)
    if len(metadata["targetGenes"]) > 1:
        msg = f"Multiple target genes for {scoreset_urn} -- look into this."
        raise ResourceAcquisitionError(msg)
    gene = metadata["targetGenes"][0]
    try:
        structured_data = ScoresetMetadata(
            urn=metadata["urn"],
            target_gene_name=gene["name"],
            target_gene_category=gene["category"],
            target_sequence=gene["targetSequence"]["sequence"],
            target_sequence_type=gene["targetSequence"]["sequenceType"],
            target_uniprot_ref=_get_uniprot_ref(metadata),
        )
    except (KeyError, ValidationError) as e:
        msg = f"Unable to extract metadata from API response for scoreset {scoreset_urn}: {e}"
        _logger.error(msg)
        raise ResourceAcquisitionError(msg) from e

    return structured_data


def _load_scoreset_records(path: Path) -> list[ScoreRow]:
    """Load scoreset records from CSV file.

    This method is intentionally identified as "private", but is refactored out for
    use during testing.
    """
    scores_data: list[ScoreRow] = []
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["score"] == "NA":
                row["score"] = None
            else:
                row["score"] = row["score"]
            scores_data.append(ScoreRow(**row))
    return scores_data


def _get_experiment_53_scores(outfile: Path, silent: bool) -> None:
    """Scores for `urn:mavedb:00000053-a-1` can be hard to acquire from the server
    on account of their considerable size. This method uses a basic workaround to fetch
    a copy hosted on a GitHub issue until we have a final resolution.
    """
    url = "https://github.com/VariantEffect/mavedb-api/files/13746791/00000053-a-1.zip"
    with tempfile.NamedTemporaryFile() as temp_file:
        path = Path(temp_file.name)
        http_download(url, path, silent)
        with zipfile.ZipFile(path, "r") as zip_ref:
            with zip_ref.open("00000053-a-1/00000053-a-1.scores.csv") as file:
                contents = file.read()
            with outfile.open("wb") as f:
                f.write(contents)


def get_scoreset_records(
    urn: str, silent: bool = True, dcd_mapping_dir: Path | None = None
) -> list[ScoreRow]:
    """Get scoreset records.

    Only hit the MaveDB API if unavailable locally. That means data must be refreshed
    manually (i.e. you'll need to delete a scoreset file yourself for this method to
    fetch a new one). This could be improved in future versions.

    :param urn: URN for scoreset
    :param silent: if true, suppress console output
    :param dcd_mapping_dir: optionally declare save location for records
    :return: Array of individual ScoreRow objects, containing information like protein
        changes and scores
    :raise ResourceAcquisitionError: if unable to fetch file
    """
    if not dcd_mapping_dir:
        dcd_mapping_dir = LOCAL_STORE_PATH
    scores_csv = dcd_mapping_dir / f"{urn}_scores.csv"
    # TODO use smarter/more flexible caching methods
    if not scores_csv.exists():
        if urn == "urn:mavedb:00000053-a-1":
            _get_experiment_53_scores(scores_csv, silent)
        else:
            url = f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores"
            try:
                http_download(url, scores_csv, silent)
            except requests.HTTPError as e:
                msg = f"HTTPError when fetching scores CSV from {url}"
                _logger.error(msg)
                raise ResourceAcquisitionError(msg) from e

    return _load_scoreset_records(scores_csv)
