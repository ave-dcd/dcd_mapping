"""Provide class definitions for commonly-used information objects."""
from decimal import Decimal
from enum import Enum
from typing import Dict, List, Literal, Optional, Union

from cool_seq_tool.schemas import AnnotationLayer, Strand, TranscriptPriority
from ga4gh.core import sha512t24u
from ga4gh.vrs._internal.models import Allele
from pydantic import BaseModel, StrictBool, StrictFloat, StrictInt, StrictStr


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
    target_uniprot_ref: Optional[UniProtRef] = None


class ScoreRow(BaseModel):
    """Row from a MAVE score result"""

    hgvs_pro: str
    hgvs_nt: str
    score: Decimal
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


class ReferenceSequence(BaseModel):
    """Base reference sequence class."""

    sequence_type: TargetSequenceType
    sequence_id: StrictStr


class ComputedReferenceSequence(ReferenceSequence):
    """Define metadata describing a computed reference sequence"""

    sequence: StrictStr


class MappedReferenceSequence(ReferenceSequence):
    """Define metadata describing a mapped, human reference sequence"""

    sequence_accessions: List[StrictStr]


class MappedOutput(BaseModel):
    """Define output format for mapped score set"""

    pre_mapped: Union[dict, List[dict]]
    post_mapped: Union[dict, List[dict]]
    mavedb_id: StrictStr
    relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"
    score: Optional[StrictFloat]


class VrsRefAlleleSeq(BaseModel):
    """Define reference sequence indicated by sequence digest
    and start and end positions in an allele
    """

    vrs_ref_allele_seq: StrictStr


class HgvsExpression(BaseModel):
    """Define class for defining an HGVS expression for a mapped
    MAVE variant
    """

    type: StrictStr = "Expression"
    syntax: StrictStr
    value: StrictStr
    syntax_version: StrictStr = None


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


class VrsObject1_x(BaseModel):  # noqa: N801
    """Define response object for VRS 1.x object"""

    mavedb_id: StrictStr
    pre_mapped_variants: Dict
    post_mapped_variants: Dict
    score: Decimal
    layer: AnnotationLayer
    relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"


class VrsMapping(BaseModel):
    """Define pre-post mapping pair structure for VRS-structured variations."""

    mavedb_id: StrictStr
    pre_mapped_protein: Optional[Union[Allele, List[Allele]]] = None
    post_mapped_protein: Optional[Union[Allele, List[Allele]]] = None
    pre_mapped_genomic: Optional[Union[Allele, List[Allele]]] = None
    post_mapped_genomic: Optional[Union[Allele, List[Allele]]] = None
    mapped_transcript: Optional[TranscriptDescription] = None
    score: Decimal
    relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"

    def serialize(self, sequence: str, start: int, end: int, sequence_id: str) -> str:
        """Get VRS 1.X Allele Digest
        :param sequence: The string given by the sequence.state attribute
        :param start: The start position given in the allele
        :param end: The end position given in the allele
        :param sequence_id: The GA4GH SQ digest
        :return A GA4GH VA digest
        """
        location_raw = f'{{"end":{{"type":"Number","value":{end}}},"sequence_id":"{sequence_id.split(".")[1]}","start":{{"type":"Number","value":{start}}},"type":"SequenceLocation"}}'
        location_serialized = sha512t24u(location_raw.encode("ascii"))
        allele_raw = f'{{"location":"{location_serialized}","state":{{"sequence":"{sequence}","type":"LiteralSequenceExpression"}},"type":"Allele"}}'
        return sha512t24u(allele_raw.encode("ascii"))

    def generate_allele_structure(self, var: Allele) -> Dict:
        """Generate VRS 1.x allele structure
        :param var: A VRS 2.0alpha allele
        :return A dictionary
        """
        sequence = "" if not var.state.sequence else var.state.sequence.root
        return {
            "type": "Allele",
            "location": {
                "id": None,
                "type": "SequenceLocation",
                "sequence_id": f"ga4gh:{var.location.sequenceReference.refgetAccession}",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {"value": var.location.start, "type": "number"},
                    "end": {"value": var.location.end, "type": "number"},
                },
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": sequence,
            },
        }

    def output_vrs_variations(self, layer: AnnotationLayer) -> VrsObject1_x:
        """Construct VRS 1.3 compatible objects from 2.0a models.

        :param layer: The Annotation Layer (genomic or protein)
        :return A VrsObject1_x object
        """
        if not self.pre_mapped_genomic and layer == AnnotationLayer.GENOMIC:
            return None
        if not self.pre_mapped_protein and layer == AnnotationLayer.PROTEIN:
            return None

        if layer == AnnotationLayer.GENOMIC:
            pre_mapped_2_0 = self.pre_mapped_genomic
            post_mapped_2_0 = self.post_mapped_genomic
        else:
            pre_mapped_2_0 = self.pre_mapped_protein
            post_mapped_2_0 = self.post_mapped_protein

        pre_mapped_variants = []
        post_mapped_variants = []

        for var in pre_mapped_2_0:
            pre_mapped = self.generate_allele_structure(var)
            pre_mapped_id = self.serialize(
                pre_mapped["state"]["sequence"],
                pre_mapped["location"]["interval"]["start"]["value"],
                pre_mapped["location"]["interval"]["end"]["value"],
                pre_mapped["location"]["sequence_id"],
            )
            pre_mapped["id"] = f"ga4gh:VA.{pre_mapped_id}"
            pre_mapped_variants.append(pre_mapped)

        for var in post_mapped_2_0:
            post_mapped = self.generate_allele_structure(var)
            post_mapped_id = self.serialize(
                post_mapped["state"]["sequence"],
                post_mapped["location"]["interval"]["start"]["value"],
                post_mapped["location"]["interval"]["end"]["value"],
                post_mapped["location"]["sequence_id"],
            )
            post_mapped["id"] = f"ga4gh:VA.{post_mapped_id}"
            post_mapped_variants.append(post_mapped)

        if len(pre_mapped_variants) > 1:
            pre_mapped_variants = {"type": "Haplotype", "members": pre_mapped_variants}
        else:
            pre_mapped_variants = pre_mapped_variants[0]

        if len(post_mapped_variants) > 1:
            post_mapped_variants = {
                "type": "Haplotype",
                "members": post_mapped_variants,
            }
        else:
            post_mapped_variants = post_mapped_variants[0]

        return VrsObject1_x(
            mavedb_id=self.mavedb_id,
            pre_mapped_variants=pre_mapped_variants,
            post_mapped_variants=post_mapped_variants,
            layer=layer,
            score=self.score,
        )
