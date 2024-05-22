"""Define Allele and Haplotype objects, as well as component parts, in conformance with
a subset of VRS 1.3 and VRSATILE.

Note that this is not an exhaustive definition of those components -- rather, the schema
generated here validates objects that are also validated by the VRS 1.3 schema,
but other valid VRS 1.3 objects might not be compatible with the schema defined here.
"""
from typing import Any, Literal

from pydantic import BaseModel, StrictInt, StrictStr


class Number(BaseModel):
    """Define VRS 1.3 Number."""

    type: Literal["number"] = "number"
    value: StrictInt


class SequenceInterval(BaseModel):
    """Define VRS 1.3 SequenceInterval."""

    type: Literal["SequenceInterval"] = "SequenceInterval"
    start: Number
    end: Number


class SequenceLocation(BaseModel):
    """Define VRS 1.3 SequenceLocation."""

    type: Literal["SequenceLocation"] = "SequenceLocation"
    sequence_id: StrictStr
    interval: SequenceInterval


class LiteralSequenceExpression(BaseModel):
    """Define VRS 1.3 LiteralSequenceExpression."""

    type: Literal["LiteralSequenceExpression"]
    sequence: StrictStr


class Allele(BaseModel):
    """Define VRS 1.3 Allele."""

    type: Literal["Allele"] = "Allele"
    location: SequenceLocation
    state: LiteralSequenceExpression


class Expression(BaseModel):
    """Define VRS 1.3 Expression."""

    type: Literal["Expression"] = "Expression"
    syntax: StrictStr
    value: StrictStr
    syntax_version: Any  # TODO ???


class VariationDescriptor(BaseModel):
    """Define VRSATILE VariationDescriptor."""

    id: StrictStr
    type: Literal["VariationDescriptor"] = "VariationDescriptor"
    # variation: Allele | Haplotype
    expresions: list[Expression]
    vrs_ref_allele_seq: StrictStr