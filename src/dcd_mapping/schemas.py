"""Provide class definitions for commonly-used information objects."""
from enum import Enum
from typing import Any, Literal

from cool_seq_tool.schemas import AnnotationLayer, Strand, TranscriptPriority
from ga4gh.core import sha512t24u
from ga4gh.vrs._internal.models import Allele, Haplotype
from pydantic import BaseModel, StrictBool, StrictInt, StrictStr

from dcd_mapping import vrs_v1_schemas


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
    # sometimes the score is "NA" or other strings, so this can't be a Decimal
    score: str
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


###################### NEW: ####################################


class MappedScore(BaseModel):
    """Provide mappings for an individual experiment score.

    This model defines the output of the VRS mapping phase of the pipeline.
    """

    accession_id: StrictStr
    annotation_layer: AnnotationLayer
    score: str
    pre_mapped: Allele | Haplotype
    post_mapped: Allele | Haplotype


class ScoreAnnotation(BaseModel):
    """Provide extra annotations on top of mappings for an individual experiment score.

    This model defines what an individual mapping instance looks like in the final JSON.
    """

    pre_mapped: vrs_v1_schemas.VariationDescriptor | vrs_v1_schemas.Haplotype
    post_mapped: vrs_v1_schemas.VariationDescriptor | vrs_v1_schemas.Haplotype
    pre_mapped_2_0: Allele | Haplotype | None = None
    post_mapped_2_0: Allele | Haplotype | None = None
    mavedb_id: StrictStr
    relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"
    score: str


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


def allele_to_v1_allele(allele: Allele) -> vrs_v1_schemas.Allele:
    """Convert VRS 2.0 allele to VRS 1.3 allele.


    :param allele: VRS 2.0a allele
    :return: equivalent VRS 1.3 allele
    """
    start = allele.location.start  # type: ignore
    end = allele.location.end  # type: ignore
    sequence_id = f"ga4gh:{allele.location.sequenceReference.refgetAccession}"  # type: ignore
    location_raw = f'{{"end":{{"type":"Number","value":{end}}},"sequence_id":"{sequence_id.split(".")[1]}","start":{{"type":"Number","value":{start}}},"type":"SequenceLocation"}}'
    location_id = sha512t24u(location_raw.encode("ascii"))
    sequence = "" if not allele.state.sequence else allele.state.sequence.root
    allele_raw = f'{{"location":"{location_id}","state":{{"sequence":"{sequence}","type":"LiteralSequenceExpression"}},"type":"Allele"}}'
    allele_id = sha512t24u(allele_raw.encode("ascii"))

    return vrs_v1_schemas.Allele(
        id=f"ga4gh:VA.{allele_id}",
        location=vrs_v1_schemas.SequenceLocation(
            id=location_id,
            sequence_id=sequence_id,
            interval=vrs_v1_schemas.SequenceInterval(
                start=vrs_v1_schemas.Number(value=start, type="number"),  # type: ignore
                end=vrs_v1_schemas.Number(value=end, type="number"),  # type: ignore
            ),
        ),
        state=vrs_v1_schemas.LiteralSequenceExpression(sequence=sequence),
    )


def allele_to_vod(allele: Allele) -> vrs_v1_schemas.VariationDescriptor:
    """Convert VRS 2.0 allele to comparable VRSATILE VariationDescriptor.

    Some allele properties aren't available in the 1.3 allele, so we have to lift them
    up to the VariationDescriptor.
    """
    allele_v1 = allele_to_v1_allele(allele)
    original_expression = allele.expressions[0]  # type: ignore
    expresions = [
        vrs_v1_schemas.Expression(
            syntax=original_expression.syntax.value,
            value=original_expression.value,
            syntax_version=None,
        )
    ]
    return vrs_v1_schemas.VariationDescriptor(
        id=allele_v1.id,
        variation=allele_v1,
        type="VariationDescriptor",
        expressions=expresions,
        vrs_ref_allele_seq=allele.extensions[0].value,  # type: ignore
        extensions=[],
    )


def haplotype_to_haplotype_1_3(haplotype: Haplotype) -> vrs_v1_schemas.Haplotype:
    """Convert VRS 2.0 Haplotype to VRS 1.3 Haplotype.

    :param haplotype: VRS 2.0 haplotype
    :return: VRS 1.3 haplotype that contains VRSATILE variation descriptors (not alleles)
    """
    members = []
    allele: Allele
    for allele in haplotype.members:  # type: ignore
        members.append(allele_to_vod(allele))
    return vrs_v1_schemas.Haplotype(members=members)


###################### WORKING: ################################


# class VrsMapping1_3(BaseModel):
#     """Define response object for VRS 1.x object"""
#
#     mavedb_id: StrictStr
#     pre_mapped_variants: dict
#     post_mapped_variants: dict
#     score: str
#     layer: AnnotationLayer
#     relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"
#
#
# class VrsMapping(BaseModel):
#     """Define pre-post mapping pair structure for VRS-structured variations."""
#
#     mavedb_id: StrictStr
#     pre_mapped_protein: Allele | list[Allele] | None = None
#     post_mapped_protein: Allele | list[Allele] | None = None
#     pre_mapped_genomic: Allele | list[Allele] | None = None
#     post_mapped_genomic: Allele | list[Allele] | None = None
#     mapped_transcript: TranscriptDescription | None = None
#     score: str
#     relation: Literal["SO:is_homologous_to"] = "SO:is_homologous_to"
#
#     def serialize(self, sequence: str, start: int, end: int, sequence_id: str) -> str:
#         """Get VRS 1.X Allele Digest
#         :param sequence: The string given by the sequence.state attribute
#         :param start: The start position given in the allele
#         :param end: The end position given in the allele
#         :param sequence_id: The GA4GH SQ digest
#         :return A GA4GH VA digest
#         """
#         location_raw = f'{{"end":{{"type":"Number","value":{end}}},"sequence_id":"{sequence_id.split(".")[1]}","start":{{"type":"Number","value":{start}}},"type":"SequenceLocation"}}'
#         location_serialized = sha512t24u(location_raw.encode("ascii"))
#         allele_raw = f'{{"location":"{location_serialized}","state":{{"sequence":"{sequence}","type":"LiteralSequenceExpression"}},"type":"Allele"}}'
#         return sha512t24u(allele_raw.encode("ascii"))
#
#     def generate_allele_structure(self, var: Allele) -> dict:
#         """Generate VRS 1.x allele structure
#         :param var: A VRS 2.0alpha allele
#         :return A dictionary
#         """
#         sequence = "" if not var.state.sequence else var.state.sequence.root
#         return {
#             "type": "Allele",
#             "location": {
#                 "id": None,
#                 "type": "SequenceLocation",
#                 "sequence_id": f"ga4gh:{var.location.sequenceReference.refgetAccession}",
#                 "interval": {
#                     "type": "SequenceInterval",
#                     "start": {"value": var.location.start, "type": "number"},
#                     "end": {"value": var.location.end, "type": "number"},
#                 },
#             },
#             "state": {
#                 "type": "LiteralSequenceExpression",
#                 "sequence": sequence,
#             },
#         }
#
#     def output_vrs_variations(self, layer: AnnotationLayer) -> VrsMapping1_3:
#         """Construct VRS 1.3 compatible objects from 2.0a models.
#
#         :param layer: The Annotation Layer (genomic or protein)
#         :return A VrsObject1_x object
#         """
#         if not self.pre_mapped_genomic and layer == AnnotationLayer.GENOMIC:
#             msg = f"Cannot map on {AnnotationLayer.GENOMIC} for mapping with no genomic variations"
#             raise ValueError(msg)
#         if not self.pre_mapped_protein and layer == AnnotationLayer.PROTEIN:
#             msg = f"Cannot map on {AnnotationLayer.PROTEIN} for mapping with no protein variations"
#             raise ValueError(msg)
#
#         if layer == AnnotationLayer.GENOMIC:
#             pre_mapped_2_0 = self.pre_mapped_genomic
#             post_mapped_2_0 = self.post_mapped_genomic
#         else:
#             pre_mapped_2_0 = self.pre_mapped_protein
#             post_mapped_2_0 = self.post_mapped_protein
#
#         pre_mapped_variants = []
#         post_mapped_variants = []
#
#         for var in pre_mapped_2_0:
#             pre_mapped = self.generate_allele_structure(var)
#             pre_mapped_id = self.serialize(
#                 pre_mapped["state"]["sequence"],
#                 pre_mapped["location"]["interval"]["start"]["value"],
#                 pre_mapped["location"]["interval"]["end"]["value"],
#                 pre_mapped["location"]["sequence_id"],
#             )
#             pre_mapped["id"] = f"ga4gh:VA.{pre_mapped_id}"
#             pre_mapped_variants.append(pre_mapped)
#
#         for var in post_mapped_2_0:
#             post_mapped = self.generate_allele_structure(var)
#             post_mapped_id = self.serialize(
#                 post_mapped["state"]["sequence"],
#                 post_mapped["location"]["interval"]["start"]["value"],
#                 post_mapped["location"]["interval"]["end"]["value"],
#                 post_mapped["location"]["sequence_id"],
#             )
#             post_mapped["id"] = f"ga4gh:VA.{post_mapped_id}"
#             post_mapped_variants.append(post_mapped)
#
#         if len(pre_mapped_variants) > 1:
#             pre_mapped_variants = {"type": "Haplotype", "members": pre_mapped_variants}
#         else:
#             pre_mapped_variants = pre_mapped_variants[0]
#
#         if len(post_mapped_variants) > 1:
#             post_mapped_variants = {
#                 "type": "Haplotype",
#                 "members": post_mapped_variants,
#             }
#         else:
#             post_mapped_variants = post_mapped_variants[0]
#
#         return VrsMapping1_3(
#             mavedb_id=self.mavedb_id,
#             pre_mapped_variants=pre_mapped_variants,
#             post_mapped_variants=post_mapped_variants,
#             layer=layer,
#             score=self.score,
#         )
