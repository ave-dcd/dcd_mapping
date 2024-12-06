"""Provide core MaveDB mapping methods."""

import logging
import os
import subprocess
from importlib.metadata import version
from pathlib import Path

import click
from cool_seq_tool.schemas import AnnotationLayer

from dcd_mapping.align import AlignmentError, BlatNotFoundError, align
from dcd_mapping.annotate import annotate, save_mapped_output_json
from dcd_mapping.lookup import check_gene_normalizer, check_seqrepo, check_uta
from dcd_mapping.mavedb_data import get_scoreset_metadata, get_scoreset_records
from dcd_mapping.resource_utils import ResourceAcquisitionError
from dcd_mapping.schemas import (
    ScoreRow,
    ScoresetMetadata,
)
from dcd_mapping.transcripts import TxSelectError, select_transcript
from dcd_mapping.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


def _emit_info(msg: str, silent: bool, log_level: int = logging.INFO) -> None:
    if not silent:
        click.echo(msg)
    if log_level == logging.INFO:
        _logger.info(msg)
    elif log_level == logging.ERROR:
        _logger.error(msg)
    else:
        msg = f"Unexpected log level requested: {log_level}"
        raise ValueError(msg)


async def _check_data_prereqs(silent: bool) -> None:
    """Non-exhaustive check that data prereqs are properly configured and available."""
    _emit_info("Checking data prereqs....", silent)
    success = True
    try:
        await check_uta()
    except Exception:
        success = False
        cst_version = version("cool-seq-tool")
        _emit_info(
            f"* UTA appears to be unavailable. Check the logs for more information. For troubleshooting, we recommend checking the UTA readme (https://github.com/biocommons/uta?tab=readme-ov-file#installing-uta-locally) and the Cool-Seq-Tool installation instructions (https://coolseqtool.readthedocs.io/{cst_version}/install.html#set-up-uta). Remember that the UTA connection is configurable via a libpq URI provided under the environment variable UTA_DB_URL (see Cool-Seq-Tool docs: https://coolseqtool.readthedocs.io/{cst_version}/usage.html#environment-configuration) -- otherwise, by default it attempts a connection to `postgresql://uta_admin:uta@localhost:5433/uta/uta_20210129b`.",
            silent,
            logging.ERROR,
        )
    try:
        check_seqrepo()
    except Exception:
        success = False
        _emit_info(
            "* SeqRepo appears inaccessible or unusable. Check the logs for more information. Ensure that a local SeqRepo snapshot has been downloaded (it should've taken a while -- see https://github.com/biocommons/biocommons.seqrepo?tab=readme-ov-file#requirements), that it's located either at `/usr/local/share/seqrepo/latest` or at the location designated by the `SEQREPO_ROOT_DIR` environment variable, and that it's writeable (see https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/store.rst).",
            silent,
            logging.ERROR,
        )
    try:
        check_gene_normalizer()
    except Exception:
        success = False
        gene_norm_verison = version("gene-normalizer")
        _emit_info(
            f"* Gene Normalizer appears to be unavailable. Check the logs for more information. Note that a data snapshot needs to be acquired, or the data update routine must be called (this should've taken at least a few seconds, if not several minutes). For troubleshooting, review the Gene Normalizer installation instructions and documentation: https://gene-normalizer.readthedocs.io/{gene_norm_verison}/install.html",
            silent,
            logging.ERROR,
        )
    try:
        configured_blat_bin = os.environ.get("BLAT_BIN_PATH")
        if configured_blat_bin:
            result = subprocess.run(  # noqa: ASYNC221 S603
                configured_blat_bin,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if result.returncode == 127:
                success = False
                _emit_info(
                    f"* BLAT binary at location pointed to by BLAT_BIN_PATH env var, {configured_blat_bin}, appears to be missing. Please check that a BLAT executable is at that location.",
                    silent,
                    logging.ERROR,
                )
            elif result.returncode != 0 and result.returncode != 255:
                success = False
                _emit_info(
                    f"* BLAT binary at location pointed to by BLAT_BIN_PATH env var, {configured_blat_bin}, doesn't appear to be properly executable. Please double-check that the executable is at the proper location and has correct permissions.",
                    silent,
                    logging.ERROR,
                )
        else:
            result = subprocess.run(  # noqa: S603 ASYNC221
                "blat",  # noqa: S607
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if result.returncode == 127:
                success = False
                _emit_info(
                    "* Unable to run BLAT. The BLAT binary must be acquired separately by the user and either accessible in the $PATH or at the location pointed to by the BLAT_BIN_PATH environment variable.",
                    silent,
                    logging.ERROR,
                )
            elif result.returncode != 0 and result.returncode != 255:
                success = False
                _emit_info(
                    "* BLAT binary appear to be properly executable. Please double-check that the executable is at the proper location and has correct permissions.",
                    silent,
                    logging.ERROR,
                )
    except Exception:
        success = False
        _emit_info(
            "Encountered unexpected error while testing availability of BLAT. Check logs for more information. See README for more information about BLAT setup.",
            silent,
            logging.ERROR,
        )
    if not success:
        raise LookupError
    _emit_info("Data prereqs checks pass.", silent)


async def map_scoreset(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    output_path: Path | None = None,
    silent: bool = True,
    check_data_prereqs: bool = True,
) -> None:
    """Given information about a MAVE experiment, map to VRS and save output as JSON.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :param output_path: optional path to save output at
    :param silent: if ``True``, suppress console information output
    :param check_data_prereqs: if ``True``, check for external data availability
        before performing mapping
    """
    if check_data_prereqs:
        await _check_data_prereqs(silent)

    _emit_info(f"Performing alignment for {metadata.urn}...", silent)
    try:
        alignment_result = align(metadata, silent)
    except BlatNotFoundError as e:
        msg = "BLAT command appears missing. Ensure it is available on the $PATH or use the environment variable BLAT_BIN_PATH to point to it. See instructions in the README prerequisites section for more."
        _emit_info(msg, silent, logging.ERROR)
        raise e
    except AlignmentError as e:
        _emit_info(
            f"Alignment failed for scoreset  {metadata.urn} {e}", silent, logging.ERROR
        )
        raise e
    _emit_info("Alignment complete.", silent)

    _emit_info("Selecting reference sequence...", silent)

    # We only need to select a transcript if a score set has only protein variants.
    # If not, this step can be skipped.
    transcript = None
    preferred_layer = None
    if "urn:mavedb:00000097" in records[0].accession:
        preferred_layer = AnnotationLayer.PROTEIN
    else:
        preferred_layer = AnnotationLayer.PROTEIN
        for row in records:
            if row.hgvs_nt != "NA":
                preferred_layer = AnnotationLayer.GENOMIC
                break

    if preferred_layer == AnnotationLayer.PROTEIN:
        try:
            transcript = await select_transcript(metadata, records, alignment_result)
        except TxSelectError as e:
            _emit_info(
                f"Transcript selection failed for scoreset {metadata.urn}",
                silent,
                logging.ERROR,
            )
            raise e
        _emit_info("Reference selection complete.", silent)

    _emit_info("Mapping to VRS...", silent)
    try:
        vrs_results = vrs_map(metadata, alignment_result, records, transcript, silent)
    except VrsMapError as e:
        _emit_info(
            f"VRS mapping failed for scoreset {metadata.urn}", silent, logging.ERROR
        )
        raise e
    if vrs_results is None:
        _emit_info(f"No mapping available for {metadata.urn}", silent)
        return
    _emit_info("VRS mapping complete.", silent)

    _emit_info("Annotating metadata and saving to file...", silent)
    vrs_results = annotate(vrs_results, transcript, metadata, alignment_result)
    final_output = save_mapped_output_json(
        metadata.urn,
        vrs_results,
        alignment_result,
        transcript,
        output_path,
    )
    _emit_info(f"Annotated scores saved to: {final_output}.", silent)


async def map_scoreset_urn(
    urn: str,
    output_path: Path | None = None,
    silent: bool = True,
    check_data_prereqs: bool = True,
) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :param output_path: optional path to save output at
    :param silent: if True, suppress console information output
    :param check_data_prereqs: if ``True``, check for external data availability
    before performing mapping
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, silent)
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        raise e
    await map_scoreset(metadata, records, output_path, silent, check_data_prereqs)
