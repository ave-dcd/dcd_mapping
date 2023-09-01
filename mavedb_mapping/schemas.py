"""Provide class definitions for commonly-used information objects."""
from enum import Enum

from pydantic import BaseModel


class TargetGeneCategory(str, Enum):
    """Define target gene category options"""

    PROTEIN_CODING = "Protein coding"


class ScoresetMetadata(BaseModel):
    """Salient metadata for a scoreset."""

    urn: str
    target_gene_name: str
    target_gene_category: TargetGeneCategory
    target_sequence: str
    target_sequence_type: str


class ScoreRow(BaseModel):
    """Individual result in a scoreset."""

    hgvs_pro: str
    hgvs_nt: str
    scores: str
    accession: str
