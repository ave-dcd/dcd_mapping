"""Provide class definitions for commonly-used information objects."""

from enum import Enum
from typing import Any, Literal

from cool_seq_tool.schemas import AnnotationLayer, Strand, TranscriptPriority
from ga4gh.vrs.models import Allele, CisPhasedBlock
from pydantic import BaseModel, StrictBool, StrictInt, StrictStr


class TargetSequenceType(str, Enum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class TargetType(str, Enum):
    """Define target gene types."""

    PROTEIN_CODING = "Protein coding"
    REGULATORY = "Regulatory"
    OTHER_NC = "Other noncoding"


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
    target_uniprot_ref: UniProtRef | None = None


class ScoreRow(BaseModel):
    """Row from a MAVE score result"""

    hgvs_pro: str
    hgvs_nt: str
    score: str | None
    accession: str


class SequenceRange(BaseModel):
    """Define range over a sequence. Useful for expressing alignment query and hit results."""

    start: int
    end: int


class GeneLocation(BaseModel):
    """Gene location info, gathered from normalizer result. Likely to be incomplete."""

    chromosome: str | None = None
    start: int | None = None
    end: int | None = None


class ReferenceSequence(BaseModel):
    """Base reference sequence class."""

    sequence_type: TargetSequenceType
    sequence_id: StrictStr


class ComputedReferenceSequence(ReferenceSequence):
    """Define metadata describing a computed reference sequence"""

    sequence: StrictStr


class MappedReferenceSequence(ReferenceSequence):
    """Define metadata describing a mapped, human reference sequence"""

    sequence_accessions: list[StrictStr]


class AlignmentResult(BaseModel):
    """Define BLAT alignment output."""

    chrom: str
    strand: Strand
    coverage: float
    ident_pct: float
    query_range: SequenceRange
    query_subranges: list[SequenceRange]
    hit_range: SequenceRange
    hit_subranges: list[SequenceRange]


class TranscriptDescription(BaseModel):
    """Structured transcript description.

    Provides less information than the MANE results, but should convey what we need.
    """

    refseq_nuc: str
    refseq_prot: str
    transcript_priority: TranscriptPriority


class ManeDescription(TranscriptDescription):
    """Structured MANE data retrieval result."""

    ncbi_gene_id: str
    ensembl_gene_id: str
    hgnc_gene_id: str
    symbol: str
    name: str
    ensembl_nuc: str
    ensembl_prot: str
    grch38_chr: str
    chr_start: int
    chr_end: int
    chr_strand: str


class TxSelectResult(BaseModel):
    """Define response object from transcript selection process."""

    nm: str | None = None
    np: str
    start: StrictInt
    is_full_match: StrictBool
    transcript_mode: TranscriptPriority | None = None
    sequence: str


class MappedScore(BaseModel):
    """Provide mappings for an individual experiment score.

    This model defines the output of the VRS mapping phase of the pipeline.
    """

    accession_id: StrictStr
    annotation_layer: AnnotationLayer
    score: str | None
    pre_mapped: Allele | CisPhasedBlock
    post_mapped: Allele | CisPhasedBlock


class ScoreAnnotation(BaseModel):
    """Provide extra annotations on top of mappings for an individual experiment score.

    This model defines what an individual mapping instance looks like in the final JSON.
    """

    pre_mapped: CisPhasedBlock | Allele
    post_mapped: CisPhasedBlock | Allele | None
    mavedb_id: StrictStr
    relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"
    score: float | None


class ScoreAnnotationWithLayer(ScoreAnnotation):
    """Couple annotations with an easily-computable definition of the annotation layer
    from which they originate.

    Used for filtering individual annotations just before saving the final JSON product.
    """

    annotation_layer: AnnotationLayer


class ScoresetMapping(BaseModel):
    """Provide all mapped scores for a scoreset."""

    metadata: Any  # TODO get exact MaveDB metadata structure?
    computed_reference_sequence: ComputedReferenceSequence
    mapped_reference_sequence: MappedReferenceSequence
    mapped_scores: list[ScoreAnnotation]
