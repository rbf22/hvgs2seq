from pydantic import BaseModel
from typing import List, Literal, Optional, Dict

from .parse import VariantNorm


class EditPlan(BaseModel):
    """
    Defines the plan for applying variants to haplotypes, including policies for handling overlaps.
    """
    haplotypes: List[List[VariantNorm]]  # One list of variants per haplotype
    policy: Literal["order_by_pos", "reject_overlaps", "left_shift"]
    warnings: List[str]


class SequenceBundle(BaseModel):
    """
    A comprehensive bundle containing reference sequences, edited sequences, annotations, and provenance.
    """
    cdna_ref: str
    cdna_edited: List[str]  # One edited sequence per haplotype
    mrna_ref: str
    mrna_edited: List[str]
    protein_ref: Optional[str] = None
    protein_edited: List[Optional[str]] = []
    annotations: Dict  # For consequences, frameshifts, stops, etc.
    provenance: Dict  # For versions, transcript info, variants applied, etc.