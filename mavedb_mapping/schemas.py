"""Provide class definitions for commonly-used information objects."""
from enum import Enum

from pydantic import BaseModel


class TargetGeneCategory(str, Enum):
    """Define target gene category options. Add more definitions as needed."""

    PROTEIN_CODING = "Protein coding"


class TargetSequenceType(str, Enum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"


class ScoresetMetadata(BaseModel):
    """Salient metadata for a scoreset."""

    urn: str
    target_gene_name: str
    target_gene_category: TargetGeneCategory
    target_sequence: str
    target_sequence_type: TargetSequenceType


class ScoreRow(BaseModel):
    """Individual result in a scoreset."""

    hgvs_pro: str
    hgvs_nt: str
    score: str
    accession: str
