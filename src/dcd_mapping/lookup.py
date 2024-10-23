"""Handle API lookups to external (non-MaveDB) services.

Data sources/handlers include:

* `CoolSeqTool <https://github.com/GenomicMedLab/cool-seq-tool/>`_
* `Gene Normalizer <https://github.com/cancervariants/gene-normalization>`_
* the `VRS-Python Translator tool <https://github.com/ga4gh/vrs-python>`_
* the UniProt web API
"""

import logging
import os
import re

import polars as pl
import requests
from biocommons.seqrepo import SeqRepo
from biocommons.seqrepo.seqaliasdb.seqaliasdb import sqlite3
from cool_seq_tool.app import (
    CoolSeqTool,
)
from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.schemas import TranscriptPriority
from ga4gh.core.domain_models import Gene
from ga4gh.core.entity_models import Extension
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.models import (
    Allele,
    LiteralSequenceExpression,
    SequenceLocation,
)
from gene.database import create_db
from gene.query import QueryHandler
from gene.schemas import MatchType, SourceName

from dcd_mapping.schemas import GeneLocation, ManeDescription, ScoresetMetadata

__all__ = [
    "CoolSeqToolBuilder",
    "get_seqrepo",
    "GeneNormalizerBuilder",
    "get_protein_accession",
    "get_transcripts",
    "get_gene_symbol",
    "get_gene_location",
    "get_chromosome_identifier",
    "get_ucsc_chromosome_name",
    "get_chromosome_identifier_from_vrs_id",
    "get_sequence",
    "translate_hgvs_to_vrs",
    "get_mane_transcripts",
    "get_uniprot_sequence",
]
_logger = logging.getLogger(__name__)

# ---------------------------------- Global ---------------------------------- #


class DataLookupError(Exception):
    """Raise for misc. issues related to resource acquisition/lookup."""


class CoolSeqToolBuilder:
    """Singleton constructor for ``cool-seq-tool`` instance."""

    def __new__(cls) -> CoolSeqTool:
        """Provide ``CoolSeqTool`` instance. Construct it if unavailable.

        :return: singleton instance of CoolSeqTool
        """
        if not hasattr(cls, "instance"):
            root_dir = os.environ.get(
                "SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest"
            )
            sr = SeqRepo(root_dir, writeable=True)
            cls.instance = CoolSeqTool(sr=sr)

        return cls.instance


def get_seqrepo() -> SeqRepoAccess:
    """Retrieve SeqRepo access instance."""
    cst = CoolSeqToolBuilder()
    return cst.seqrepo_access


class GeneNormalizerBuilder:
    """Singleton constructor for Gene Normalizer instance."""

    def __new__(cls) -> QueryHandler:
        """Provide Gene Normalizer instance. Construct it if unavailable.

        :return: singleton instance of ``QueryHandler`` for Gene Normalizer
        """
        if not hasattr(cls, "instance"):
            db = create_db()
            q = QueryHandler(db)
            cls.instance = q
        return cls.instance


class TranslatorBuilder:
    """Singleton constructor for VRS Translator instance."""

    def __new__(cls, data_proxy: SeqRepoDataProxy) -> AlleleTranslator:
        """Provide translator instance. Constructs it if unavailable. Use a new
        ``data_proxy`` instance that contains a given score row's sequence/ID.

        :return: singleton instance of ``AlleleTranslator``
        """
        if not hasattr(cls, "instance"):
            tr = AlleleTranslator(data_proxy)
            cls.instance = tr
        else:
            cls.instance.data_proxy = data_proxy
        return cls.instance


# ----------------------------------- UTA ----------------------------------- #


async def check_uta() -> None:
    """Check that UTA connection appears to be working.

    :raise LookupError: if schema check fails. Realistically, if UTA isn't working right,
        another error will probably get raised first.
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"select * from {uta.schema}.meta"  # noqa: S608
    result = await uta.execute_query(query)
    if not result:
        raise DataLookupError


async def get_protein_accession(transcript: str) -> str | None:
    """Retrieve protein accession for a transcript.

    :param transcript: transcript accession, e.g. ``"NM_002529.3"``
    :return: protein accession if successful
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"""
    SELECT pro_ac FROM {uta.schema}.associated_accessions
    WHERE tx_ac = '{transcript}'
    """  # noqa: S608
    result = await uta.execute_query(query)
    if result:
        return result[0]["pro_ac"]
    return None


async def check_gene_uta(gene: str) -> bool:
    """Check if an HGNC gene symbol exists in UTA

    :param gene: An HGNC gene symbol
    :return: True if the gene symbol exists in UTA, False if not
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"""
    SELECT DISTINCT HGNC FROM {uta.schema}.tx_exon_aln_v
    WHERE hgnc = '{gene}'
    """  # noqa: S608
    result = await uta.execute_query(query)
    return bool(result)


async def get_transcripts(
    gene_symbol: str, chromosome_ac: str, start: int, end: int
) -> list[str]:
    """Get transcript accessions matching given parameters (excluding non-coding RNA).

    TODO: may be able to successfully query with only one of gene symbol/chromosome ac.
    In initial testing, gene symbol doesn't seem to be a meaningful filter, but should
    get further confirmation.

    :param gene_symbol: HGNC-given gene symbol (usually, but not always, equivalent to
        symbols available in other nomenclatures.)
    :param chromosome: chromosome accession (e.g. ``"NC_000007.13"``)
    :param start: starting position
    :param end: ending position
    :return: candidate transcript accessions
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"""
    SELECT tx_ac
    FROM {uta.schema}.tx_exon_aln_v
    WHERE hgnc = '{gene_symbol}'
      AND ({start} BETWEEN alt_start_i AND alt_end_i OR {end} BETWEEN alt_start_i AND
      alt_end_i OR (alt_start_i >= {start} AND alt_end_i <= {end}))
      AND alt_ac = '{chromosome_ac}'
      AND tx_ac NOT LIKE 'NR_%';
    """  # noqa: S608
    result = await uta.execute_query(query)
    return [row["tx_ac"] for row in result]


# ------------------------------ Gene Normalizer ------------------------------ #


def check_gene_normalizer() -> None:
    q = GeneNormalizerBuilder()
    if (not q.db.check_schema_initialized()) or not (q.db.check_tables_populated()):
        raise DataLookupError
    if q.normalize("BRAF").match_type == MatchType.NO_MATCH:
        raise DataLookupError


async def _get_hgnc_symbol(term: str) -> str | None:
    """Fetch HGNC symbol from gene term.

    :param term: gene referent
    :return: gene symbol if available
    """
    q = GeneNormalizerBuilder()
    result = q.normalize_unmerged(term)
    hgnc = result.source_matches.get(SourceName.HGNC)
    if hgnc and len(hgnc.records) > 0:
        # probably fine to just use first match
        symbol = hgnc.records[0].symbol
        uta_status = await check_gene_uta(symbol)
        if uta_status is True:
            return symbol
        # Check if previous symbol occurs in UTA
        for symbol in hgnc.records[0].previous_symbols:
            uta_status = await check_gene_uta(symbol)
            if uta_status is True:
                return symbol
    return None


async def get_gene_symbol(metadata: ScoresetMetadata) -> str | None:
    """Acquire HGNC gene symbol given provided metadata from scoreset.

    Right now, we use two sources for normalizing:
    1. UniProt ID, if available
    2. Target name: specifically, we try the first word in the name (this could
    cause some problems and we should double-check it)

    :param ScoresetMetadata: data given by MaveDB API
    :return: gene symbol if available
    """
    if metadata.target_uniprot_ref:
        result = await _get_hgnc_symbol(metadata.target_uniprot_ref.id)
        if result:
            return result

    # try taking the first word in the target name
    if metadata.target_gene_name:
        parsed_name = ""
        if "_" in metadata.target_gene_name:
            parsed_name = metadata.target_gene_name.split("_")
        elif " " in metadata.target_gene_name:
            parsed_name = metadata.target_gene_name.split(" ")
        else:
            parsed_name = [metadata.target_gene_name]
        for word in parsed_name:
            gene = await _get_hgnc_symbol(word)
            if gene:
                return gene
    return None


def _normalize_gene(term: str) -> Gene | None:
    """Fetch normalizer response for gene term.

    :param term: gene name or referent to normalize
    :return: normalized Gene if successful
    """
    q = GeneNormalizerBuilder()
    response = q.normalize(term)
    if response.match_type > 0:
        return response.gene
    return None


def get_normalized_gene_response(
    metadata: ScoresetMetadata,
) -> Gene | None:
    """Fetch best normalized concept given available scoreset metadata.

    :param metadata: salient scoreset metadata items
    :return: Normalized gene if available
    """
    if metadata.target_uniprot_ref:
        gene_descriptor = _normalize_gene(metadata.target_uniprot_ref.id)
        if gene_descriptor:
            return gene_descriptor

    # try taking the first word in the target name
    if metadata.target_gene_name and "Minigene" not in metadata.target_gene_name:
        parsed_name = ""
        if "_" in metadata.target_gene_name:
            parsed_name = metadata.target_gene_name.split("_")
        elif " " in metadata.target_gene_name:
            parsed_name = metadata.target_gene_name.split(" ")
        else:
            parsed_name = [metadata.target_gene_name]
        for word in parsed_name:
            word = extract_potential_gene_symbol(word)
            gene_descriptor = _normalize_gene(word)
            if gene_descriptor:
                return gene_descriptor

    return None


def extract_potential_gene_symbol(word: str) -> str:
    """Remove paranetheses from a potential gene symbol

    :param word: A potential gene symbol
    :return A potential gene symbol with parantheses extracted
    """
    return re.sub(r"[(),]", "", word)


def _get_genomic_interval(
    extensions: list[Extension], src_name: str
) -> GeneLocation | None:
    """Extract start/end coords from extension list. Extensions in normalized genes
    can be of several different types, but we only want SequenceLocation data.

    :param extensions: extensions given in a descriptor
    :return: genomic interval if available
    """
    locations = [ext for ext in extensions if f"{src_name}_locations" in ext.name]
    if locations and len(locations[0].value) > 0:
        location_values = [
            v for v in locations[0].value if v["type"] == "SequenceLocation"
        ]
        if location_values:
            return GeneLocation(
                start=location_values[0]["start"],
                end=location_values[0]["end"],
                chromosome=get_chromosome_identifier_from_vrs_id(
                    f"ga4gh:{location_values[0]['sequenceReference']['refgetAccession']}"
                ),
            )
    return None


def get_gene_location(metadata: ScoresetMetadata) -> GeneLocation | None:
    """Acquire gene location data from gene normalizer using metadata provided by
    scoreset.

    As with ``get_gene_symbol()``, we try to normalize from the following:
    1. UniProt ID, if available
    2. Target name: specifically, we try the first word in the name (this could
    cause some problems and we should double-check it)

    :param metadata: data given by MaveDB API
    :return: gene location data if available
    """
    gene_descriptor = get_normalized_gene_response(metadata)
    if not gene_descriptor or not gene_descriptor.extensions:
        return None

    for src_name in ("ensembl", "ncbi"):
        loc = _get_genomic_interval(gene_descriptor.extensions, src_name)
        if loc:
            return loc

    return None


# --------------------------------- SeqRepo --------------------------------- #


def check_seqrepo() -> None:
    sr = get_seqrepo()
    if not sr.sr["NC_000001.11"][780000:780020]:
        raise DataLookupError
    try:
        conn = sr.sr.aliases._db  # noqa: SLF001
        cursor = conn.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS test_table (id INTEGER PRIMARY KEY)")
        cursor.execute("INSERT INTO test_table (id) VALUES (1)")
        conn.commit()
        cursor.execute("DELETE FROM test_table WHERE id = 1")
        cursor.execute("DROP TABLE test_table")
        conn.commit()
        conn.close()
    except sqlite3.Error as e:
        _logger.error("SeqRepo sequences DB isn't writeable.")
        raise DataLookupError from e
    # delete the CST instance so that future calls are made against an open SeqAlias
    # connection, since we close it above
    del CoolSeqToolBuilder.instance


def get_chromosome_identifier(chromosome: str) -> str:
    """Get latest NC_ accession identifier given a chromosome name.

    :param chromosome: chromosome name, e.g. ``"8"``, ``"X"``
    :return: latest ID if available
    :raise KeyError: if unable to retrieve identifier
    """
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"
    sr = get_seqrepo()
    acs = []
    for assembly in ["GRCh38", "GRCh37"]:
        tmp_acs, _ = sr.translate_identifier(
            f"{assembly}:{chromosome}", target_namespaces="refseq"
        )
        for ac in tmp_acs:
            acs.append(ac.split("refseq:")[-1])  # noqa: PERF401
    if not acs:
        raise KeyError

    # make sure e.g. version .10 > version .9
    sorted_results = sorted(acs, key=lambda i: int(i.split(".")[-1]))
    return sorted_results[-1]


def get_ucsc_chromosome_name(chromosome: str) -> str:
    """Get UCSC/GENCODE-style chromosome name, eg ``"chr1"`` instead of ``"1"`` or
    ``"NC_000001.11"``.

    :param chromosome: chromosome name/identifier
    :return: UCSC/GENCODE-style chromosome name
    :raise KeyError: if unable to find matching name
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(chromosome, "GRCh38")
    if not result:
        raise KeyError

    sorted_results = sorted([r for r in result if "chr" in r])
    try:
        return sorted_results[-1].split(":")[1]
    except IndexError as e:
        raise KeyError from e


def get_chromosome_identifier_from_vrs_id(sequence_id: str) -> str | None:
    """Get NC_ identifier given a VRS sequence ID.

    :param sequence_id: identifier a la ``ga4gh:SQ.XXXXXX``
    :return: NC_ chromosome ID
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(sequence_id, "refseq")
    if not result:
        raise KeyError

    sorted_results = sorted(result)
    return sorted_results[-1]


def get_vrs_id_from_identifier(sequence_id: str) -> str | None:
    """Get GA4GH SQ identifier given an NP_ sequence id:
    :param: GA4GH SQ digest
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(sequence_id, "ga4gh")
    if not result:
        raise KeyError

    sorted_results = sorted(result)
    return sorted_results[-1]


def get_sequence(
    sequence_id: str,
    start: int | None = None,
    end: int | None = None,
) -> str:
    """Get reference sequence given a sequence identifier.

    :param sequence_id: sequence identifier, e.g. ``"NP_938033.1"``
    :return: sequence
    :raise KeyError: if lookup fails
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    try:
        sequence = sr.get_sequence(sequence_id, start, end)
    except (KeyError, ValueError) as e:
        _logger.error("Unable to acquire sequence for ID: %s", sequence_id)
        raise KeyError from e
    if sequence is None:
        _logger.error("Unable to acquire sequence for ID: %s", sequence_id)
        raise KeyError
    return sequence


# -------------------------------- VRS-Python -------------------------------- #


def translate_hgvs_to_vrs(hgvs: str) -> Allele | None:
    """Convert HGVS variation description to VRS object.

    :param hgvs: MAVE-HGVS variation string
    :return: Corresponding VRS allele as a Pydantic class
    """
    # coerce tmp HGVS string into formally correct term
    if hgvs.startswith("NC_") and ":c." in hgvs:
        hgvs = hgvs.replace(":c.", ":g.")

    tr = TranslatorBuilder(get_seqrepo())
    try:
        allele: Allele = tr.translate_from(hgvs, "hgvs", do_normalize=False)
    except ValueError as warn:
        _logger.warning(warn)
        return None

    if (
        not isinstance(allele.location, SequenceLocation)
        or not isinstance(allele.location.start, int)
        or not isinstance(allele.location.end, int)
        or not isinstance(allele.state, LiteralSequenceExpression)
    ):
        raise ValueError
    return allele


# ----------------------------------- MANE ----------------------------------- #


def get_mane_transcripts(transcripts: list[str]) -> list[ManeDescription]:
    """Get corresponding MANE data for transcripts. Results given in order of
    transcript preference.

    :param transcripts: candidate transcripts list
    :return: complete MANE descriptions
    """

    def _sort_mane_result(description: ManeDescription) -> int:
        if description.transcript_priority == TranscriptPriority.MANE_SELECT:
            return 2
        if description.transcript_priority == TranscriptPriority.MANE_PLUS_CLINICAL:
            return 1
        # should be impossible
        _logger.warning(
            "Unrecognized transcript priority value %s for transcript description of %s",
            description.transcript_priority,
            description.refseq_nuc,
        )
        return 0

    mane_df = CoolSeqToolBuilder().mane_transcript_mappings.df
    mane_results = mane_df.filter(pl.col("RefSeq_nuc").is_in(transcripts))
    mane_data = [
        ManeDescription(
            ncbi_gene_id=row["#NCBI_GeneID"],
            ensembl_gene_id=row["Ensembl_Gene"],
            hgnc_gene_id=row["HGNC_ID"],
            symbol=row["symbol"],
            name=row["name"],
            refseq_nuc=row["RefSeq_nuc"],
            refseq_prot=row["RefSeq_prot"],
            ensembl_nuc=row["Ensembl_nuc"],
            ensembl_prot=row["Ensembl_prot"],
            transcript_priority=TranscriptPriority(
                "_".join(row["MANE_status"].lower().split())
            ),
            grch38_chr=row["GRCh38_chr"],
            chr_start=row["chr_start"],
            chr_end=row["chr_end"],
            chr_strand=row["chr_strand"],
        )
        for row in mane_results.rows(named=True)
    ]
    mane_data.sort(key=_sort_mane_result)
    return mane_data


# ---------------------------------- Misc. ---------------------------------- #


def get_uniprot_sequence(uniprot_id: str) -> str | None:
    """Get sequence directly from UniProt.

    :param uniprot_id: ID provided with target info
    :return: transcript accession if successful
    :raise HTTPError: if response comes with an HTTP error code
    """
    url = f"https://www.ebi.ac.uk/proteins/api/proteins?accession={uniprot_id.split(':')[1]}&format=json"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    json = response.json()
    return json[0]["sequence"]["sequence"]
