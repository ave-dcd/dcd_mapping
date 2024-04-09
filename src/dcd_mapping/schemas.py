"""Provide class definitions for commonly-used information objects."""
from enum import StrEnum
from typing import Dict, List, Optional, Tuple, Union

from cool_seq_tool.schemas import AnnotationLayer, Strand, TranscriptPriority
from ga4gh.vrs._internal.models import Allele
from pydantic import BaseModel, StrictBool, StrictFloat, StrictInt, StrictStr


class TargetSequenceType(StrEnum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class ReferenceGenome(StrEnum):
    """Define known reference genome names."""

    HG38 = "hg38"
    HG19 = "hg19"
    HG16 = "hg16"


class TargetType(StrEnum):
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
    target_reference_genome: ReferenceGenome
    target_uniprot_ref: Optional[UniProtRef] = None


class ScoreRow(BaseModel):
    """Row from a MAVE score result"""

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
    """Define BLAT alignment output."""

    chrom: str
    strand: Strand
    coverage: float
    ident_pct: float
    query_range: SequenceRange
    query_subranges: List[SequenceRange]
    hit_range: SequenceRange
    hit_subranges: List[SequenceRange]


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

    nm: Optional[str] = None
    np: str
    start: StrictInt
    is_full_match: StrictBool
    transcript_mode: Optional[TranscriptPriority] = None
    sequence: str


class VrsMapping(BaseModel):
    """Define pre-post mapping pair structure for VRS-structured variations."""

    mavedb_id: StrictStr
    pre_mapped_protein: Optional[Union[Allele, List[Allele]]] = None
    post_mapped_protein: Optional[Union[Allele, List[Allele]]] = None
    pre_mapped_genomic: Optional[Union[Allele, List[Allele]]] = None
    post_mapped_genomic: Optional[Union[Allele, List[Allele]]] = None
    mapped_transcript: Optional[TranscriptDescription] = None
    score: Union[StrictFloat, str]
    # relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"

    def output_vrs_variations(self, layer: AnnotationLayer) -> Tuple[Dict, Dict]:
        """Construct VRS 1.3 compatible objects from 2.0a models."""
        if layer == AnnotationLayer.GENOMIC:
            pre_mapped_2_0 = self.pre_mapped_genomic
            post_mapped_2_0 = self.post_mapped_genomic
        else:  # protein coding
            pre_mapped_2_0 = self.pre_mapped_protein
            post_mapped_2_0 = self.post_mapped_protein

        # TODO do we need to think about haplotype?
        if not isinstance(pre_mapped_2_0, Allele) or not isinstance(
            post_mapped_2_0, Allele
        ):
            raise NotImplementedError

        pre_mapped = {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequence_id": f"ga4gh:{pre_mapped_2_0.location.sequenceReference.refgetAccession}",
                "start": {"value": pre_mapped_2_0.location.start, "type": "number"},
                "end": {"value": pre_mapped_2_0.location.end, "type": "number"},
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": pre_mapped_2_0.state.sequence,
            },
        }
        post_mapped = {
            "type": "Allele",
            "location": {
                "type": "SequenceLocation",
                "sequence_id": f"ga4gh:{post_mapped_2_0.location.sequenceReference.refgetAccession}",
                "start": {"value": post_mapped_2_0.location.start, "type": "number"},
                "end": {"value": post_mapped_2_0.location.end, "type": "number"},
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": post_mapped_2_0.state.sequence,
            },
        }

        # run ga4gh identify

        return (pre_mapped, post_mapped)


class VrsMappingResult(BaseModel):
    """Define response object from VRS mappings method.

    Might not be necessary (should just be list of VrsMappings?)
    """

    variations: List[VrsMapping]
