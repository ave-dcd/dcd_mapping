"""Handle API lookups to external (non-MaveDB) services.

This module should contain methods that we don't want to think about caching.
"""
import logging
from typing import List, Optional

import requests
from cool_seq_tool import CoolSeqTool
from ga4gh.vrsatile.pydantic.vrsatile_models import Extension, GeneDescriptor
from gene.database import create_db
from gene.query import QueryHandler
from gene.schemas import SourceName

from mavemap.schemas import GeneLocation, ManeData, ScoresetMetadata

_logger = logging.getLogger(__name__)


# ---------------------------------- Global ---------------------------------- #


class CoolSeqToolBuilder:
    """Singleton constructor for ``cool-seq-tool`` instance."""

    def __new__(cls) -> CoolSeqTool:
        """Provide ``CoolSeqTool`` instance. Construct it if unavailable.

        :return: singleton instance of CoolSeqTool
        """
        if not hasattr(cls, "instance"):
            q = QueryHandler(create_db())
            cls.instance = CoolSeqTool(gene_query_handler=q)
        return cls.instance


# ----------------------------------- UTA ----------------------------------- #


async def get_protein_accession(transcript: str) -> Optional[str]:
    """Retrieve protein accession for a transcript.

    :param transcript: transcript accession, e.g. ``"NM_002529.3"``
    :return: protein accession if successful
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"""
    SELECT pro_ac FROM {uta.schema}.associated_accessions
    WHERE tx_ac = '{transcript}'
    """
    result = await uta.execute_query(query)
    if result:
        return result[0]["pro_ac"]


async def get_transcripts(
    gene_symbol: str, chromosome_ac: str, start: int, end: int
) -> List[str]:
    """Get transcript accessions matching given parameters (excluding non-coding RNA).

    TODO: may be able to successfully query with only one of gene symbol/chromosome ac.

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
      AND ({start} BETWEEN alt_start_i AND alt_end_i OR {end} BETWEEN alt_start_i AND alt_end_i)
      AND alt_ac = '{chromosome_ac}'
      AND tx_ac NOT LIKE 'NR_%';
    """
    result = await uta.execute_query(query)
    return [row["tx_ac"] for row in result]


def get_mane_transcripts(transcripts: List[str]) -> List[ManeData]:
    """Get corresponding MANE data for transcripts.

    :param transcripts: candidate transcripts list
    :return: complete MANE descriptions
    """
    mane = CoolSeqToolBuilder().mane_transcript_mappings
    mane_transcripts = mane.get_mane_from_transcripts(transcripts)
    mane_data = []
    for result in mane_transcripts:
        mane_data.append(
            ManeData(
                ncbi_gene_id=result["#NCBI_GeneID"],
                ensembl_gene_id=result["Ensembl_Gene"],
                hgnc_gene_id=result["HGNC_ID"],
                symbol=result["symbol"],
                name=result["name"],
                refseq_nuc=result["RefSeq_nuc"],
                refseq_prot=result["RefSeq_prot"],
                ensembl_nuc=result["Ensembl_nuc"],
                ensembl_prot=result["Ensembl_prot"],
                mane_status=result["MANE_status"],
                grch38_chr=result["GRCh38_chr"],
                chr_start=result["chr_start"],
                chr_end=result["chr_end"],
                chr_strand=result["chr_strand"],
            )
        )
    return mane_data


# ------------------------------ Gene Normalizer ------------------------------ #


def _get_hgnc_symbol(term: str) -> Optional[str]:
    """Fetch HGNC symbol from gene term.

    :param term: gene referent
    :return: gene symbol if available
    """
    q = CoolSeqToolBuilder().gene_query_handler
    result = q.normalize_unmerged(term)
    hgnc = result.source_matches.get(SourceName.HGNC)
    if hgnc and len(hgnc.records) > 0:
        # probably fine to just use first match
        return hgnc.records[0].symbol
    return None


def get_gene_symbol(metadata: ScoresetMetadata) -> Optional[str]:
    """Acquire HGNC gene symbol given provided metadata from scoreset.

    Right now, we use two sources for normalizing:
    1. UniProt ID, if available
    2. Target name: specifically, we try the first word in the name (this could
    cause some problems and we should double-check it)

    :param ScoresetMetadata: data given by MaveDB API
    :return: gene symbol if available
    """
    if metadata.target_uniprot_ref:
        result = _get_hgnc_symbol(metadata.target_uniprot_ref.id)
        if result:
            return result

    # try taking the first word in the target name
    if metadata.target_gene_name:
        parsed_name = metadata.target_gene_name.split(" ")[0]
        return _get_hgnc_symbol(parsed_name)


def _normalize_gene(term: str) -> Optional[GeneDescriptor]:
    """Fetch normalizer response for gene term.

    :param term: gene name or referent to normalize
    :return: GeneDescriptor if successful
    """
    q = CoolSeqToolBuilder().gene_query_handler
    response = q.normalize(term)
    if response.match_type > 0:
        return response.gene_descriptor
    else:
        return None


def _get_normalized_gene_response(
    metadata: ScoresetMetadata,
) -> Optional[GeneDescriptor]:
    """Fetch best normalized concept given available scoreset metadata.

    :param metadata: salient scoreset metadata items
    :return: Normalized gene descriptor if available
    """
    if metadata.target_uniprot_ref:
        gene_descriptor = _normalize_gene(metadata.target_uniprot_ref.id)
        if gene_descriptor:
            return gene_descriptor

    # try taking the first word in the target name
    if metadata.target_gene_name:
        parsed_name = metadata.target_gene_name.split(" ")[0]
        gene_descriptor = _normalize_gene(parsed_name)
        if gene_descriptor:
            return gene_descriptor

    return None


def _get_genomic_interval(
    extensions: List[Extension], src_name: str
) -> Optional[GeneLocation]:
    """Extract start/end coords from extension list. Extensions in gene descriptors
    can be of many different types, but we only want SequenceLocation data.

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
                start=location_values[0]["interval"]["start"]["value"],
                end=location_values[0]["interval"]["end"]["value"],
                chromosome=get_chromosome_identifier_from_vrs_id(
                    location_values[0]["sequence_id"]
                ),
            )
    return None


def get_gene_location(metadata: ScoresetMetadata) -> Optional[GeneLocation]:
    """Acquire gene location data from gene normalizer using metadata provided by
    scoreset.

    As with ``get_gene_symbol()``, we try to normalize from the following:
    1. UniProt ID, if available
    2. Target name: specifically, we try the first word in the name (this could
    cause some problems and we should double-check it)

    :param metadata: data given by MaveDB API
    :return: gene location data if available
    """
    gene_descriptor = _get_normalized_gene_response(metadata)
    if not gene_descriptor or not gene_descriptor.extensions:
        return None

    hgnc_locations = [
        loc for loc in gene_descriptor.extensions if loc.name == "hgnc_locations"
    ]
    if hgnc_locations and len(hgnc_locations[0].value) > 0:
        return GeneLocation(chromosome=hgnc_locations[0].value[0].chr)

    for src_name in ("ensembl", "ncbi"):
        loc = _get_genomic_interval(gene_descriptor.extensions, src_name)
        if loc:
            return loc

    return None


# --------------------------------- SeqRepo --------------------------------- #
# TODO
# * these could be refactored into a single method
# * not clear if all of them are necessary
# * either way, they should all be renamed


def get_chromosome_identifier(chromosome: str) -> str:
    """Get latest NC_ identifier given a chromosome name.

    :param chromosome: prefix-free chromosome name, e.g. ``"8"``, ``"X"``
    :return: latest ID if available
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.chromosome_to_acs(chromosome)
    if not result:
        raise KeyError

    sorted_results = sorted(result)
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
    except IndexError:
        raise KeyError


def get_chromosome_identifier_from_vrs_id(sequence_id: str) -> Optional[str]:
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


def get_reference_sequence(sequence_id: str) -> str:
    """Get reference sequence given a sequnce identifier.

    :param sequence_id: sequence identifier, e.g. ``"NP_938033.1"``
    :return: sequence
    :raise KeyError: if lookup fails
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    try:
        sequence = sr.get_sequence(sequence_id)
    except (KeyError, ValueError):
        _logger.error(f"Unable to acquire sequence for ID: {sequence_id}")
        raise KeyError
    if sequence is None:
        _logger.error(f"Unable to acquire sequence for ID: {sequence_id}")
        raise KeyError
    return sequence


# ---------------------------------- Misc. ---------------------------------- #


def get_clingen_id(hgvs: str) -> Optional[str]:
    """Fetch ClinGen ID. TODO finish this.

    :param hgvs: HGVS ID todo ??
    :return: ClinGen ID if available
    :raise HTTPError: if request encounters an error
    """
    url = f"https://reg.genome.network/allele?hgvs={hgvs}"
    response = requests.get(url)
    response.raise_for_status()
    page = response.json()
    page = page["@id"]
    return page.split("/")[4]


def get_uniprot_sequence(uniprot_id: str) -> Optional[str]:
    """Get sequence directly from UniProt.

    :param uniprot_id: ID provided with target info
    :return: transcript accession if successful
    :raise HTTPError: if response comes with an HTTP error code
    """
    url = f"https://www.ebi.ac.uk/proteins/api/proteins?accession={uniprot_id.split(':')[1]}&format=json"
    response = requests.get(url)
    response.raise_for_status()
    json = response.json()
    return json[0]["sequence"]["sequence"]
