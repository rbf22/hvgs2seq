"""Data models for hgvs2seq using Pydantic."""
from typing import List, Optional, Literal, Dict, Any
from pydantic import BaseModel, Field

class TranscriptConfig(BaseModel):
    """Configuration for a single transcript."""
    transcript_id: str = Field(..., description="Transcript ID (e.g., NM_..., ENST...)")
    gene_symbol: Optional[str] = Field(None, description="Gene symbol")
    assembly: str = Field(..., description="Reference assembly (e.g., GRCh38)")
    strand: Literal[1, -1] = Field(..., description="Strand of the transcript (+1 or -1)")
    exons: List[tuple[int, int]] = Field(..., description="Genomic coordinates of exons (1-based, inclusive)")
    cds_start_c: Optional[int] = Field(None, description="cDNA index of CDS start (ATG), 1-based")
    cds_end_c: Optional[int] = Field(None, description="cDNA index of the last base of the CDS, 1-based")

class VariantIn(BaseModel):
    """Input for a single variant."""
    hgvs: str = Field(..., description="Original HGVS string")
    phase_group: Optional[int] = Field(None, description="Phase group for haplotyping (optional)")

class VariantNorm(BaseModel):
    """Normalized variant representation."""
    hgvs_c: str = Field(..., description="Normalized HGVS on transcript (c. coordinates)")
    kind: Literal["sub", "del", "ins", "delins", "dup", "inv", "mnv"] = Field(..., description="Variant type")
    c_start: int = Field(..., description="Start position in cDNA coordinates (1-based)")
    c_start_offset: int = Field(0, description="Offset from the start position (for intronic variants)")
    c_end: int = Field(..., description="End position in cDNA coordinates (1-based)")
    c_end_offset: int = Field(0, description="Offset from the end position (for intronic variants)")
    alt: str = Field(..., description="Alternate sequence (inserted/alternate bases)")
    meta: Dict[str, Any] = Field({}, description="Metadata (original string, notes, etc.)")
    phase_group: Optional[int] = Field(None, description="Phase group for haplotyping (optional)")

class EditPlan(BaseModel):
    """Plan for applying edits to a transcript."""
    haplotypes: List[List[VariantNorm]] = Field(..., description="List of variants for each haplotype")
    policy: Literal["order_by_pos", "reject_overlaps", "left_shift"] = Field(..., description="Policy for handling variant overlaps")
    warnings: List[str] = Field([], description="Warnings generated during planning")

class ProteinOutcome(BaseModel):
    """Describes the outcome on the protein sequence."""
    protein_sequence: Optional[str] = Field(None, description="The resulting protein sequence")
    hgvs_p: Optional[str] = Field(None, description="HGVS notation for the protein change")
    consequence: str = Field(..., description="Classification of the effect (e.g., missense, frameshift)")

class NMDOutcome(BaseModel):
    """Describes the NMD outcome."""
    status: Literal["likely", "escape", "not_applicable"] = Field(..., description="NMD status")
    rationale: str = Field(..., description="Explanation for the NMD status decision")

class TranscriptOutcome(BaseModel):
    """Describes the outcome for a single transcript scenario."""
    haplotype_id: int = Field(..., description="Identifier for the haplotype (0 for unphased/combined)")
    scenario_id: str = Field("baseline", description="Identifier for the splicing scenario")
    mrna_sequence: str = Field(..., description="The edited mRNA sequence")
    protein_outcome: Optional[ProteinOutcome] = Field(None, description="Outcome on the protein")
    nmd: NMDOutcome = Field(..., description="NMD analysis outcome")
    confidence: float = Field(1.0, description="Confidence score for this scenario")
    extras: Dict[str, Any] = Field({}, description="Additional annotations (e.g., uORFs, NSD risk)")

class SequenceBundle(BaseModel):
    """Final output bundle containing all sequences and annotations."""
    primary_outcomes: List[TranscriptOutcome] = Field(..., description="Primary outcomes for each haplotype")
    alternate_outcomes: List[TranscriptOutcome] = Field([], description="Alternative outcomes (e.g., from splicing)")
    provenance: Dict[str, Any] = Field(..., description="Provenance information (versions, inputs, etc.)")
    warnings: List[str] = Field([], description="List of warnings from the entire process")