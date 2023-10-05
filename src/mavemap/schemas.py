"""Provide class definitions for commonly-used information objects."""
from enum import Enum
from typing import List, Optional

from pydantic import BaseModel


class TargetGeneCategory(str, Enum):
    """Define target gene category options. Add more definitions as needed."""

    PROTEIN_CODING = "Protein coding"


class TargetSequenceType(str, Enum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class ReferenceGenome(str, Enum):
    """Define known reference genome names."""

    HG38 = "hg38"


class TargetType(str, Enum):
    """Define target gene types."""

    PROTEIN_CODING = "Protein coding"


class UniProtRef(BaseModel):
    """Store metadata associated with MaveDB UniProt reference"""

    id: str
    offset: int


class ScoresetMetadata(BaseModel):
    """Store all relevant metadata from metadata reported for scoreset by MaveDB"""

    urn: str
    target_gene_name: str
    target_gene_category: TargetType
    target_sequence: str
    target_sequence_type: TargetSequenceType
    target_reference_genome: ReferenceGenome
    target_uniprot_ref: Optional[UniProtRef] = None


class ScoreRow(BaseModel):
    """TODO"""

    hgvs_pro: str
    hgvs_nt: str
    score: str
    accession: str


class SequenceRange(BaseModel):
    """Define range over a sequence. Useful for expressing alignment query and hit results."""

    start: int
    end: int


class GeneLocation(BaseModel):
    """Gene location info, gathered from normalizer result. Likely to be incomplete."""

    chromosome: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None


class AlignmentResult(BaseModel):
    """Structured BLAT alignment output."""

    chrom: str
    strand: str
    coverage: float
    ident_pct: float
    query_range: SequenceRange
    query_subranges: List[SequenceRange]
    hit_range: SequenceRange
    hit_subranges: List[SequenceRange]


class TranscriptStatus(str, Enum):
    """Define legal MANE statuses."""

    SELECT = "MANE Select"
    PLUS_CLINICAL = "MANE Plus Clinical"
    LONGEST_COMPATIBLE = "Longest Compatible"


class ManeData(BaseModel):
    """Structured MANE data retrieval result."""

    ncbi_gene_id: str
    ensembl_gene_id: str
    hgnc_gene_id: str
    symbol: str
    name: str
    refseq_nuc: str
    refseq_prot: str
    ensembl_nuc: str
    ensembl_prot: str
    mane_status: TranscriptStatus
    grch38_chr: str
    chr_start: int
    chr_end: int
    chr_strand: str
