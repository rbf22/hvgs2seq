"""Data models for hgvs2seq using Pydantic."""
from enum import Enum
from typing import List, Optional, Literal, Dict, Any, Union, Tuple
from pydantic import BaseModel, Field, field_validator, model_validator
from pathlib import Path
import os
from dataclasses import dataclass


class VariantType(str, Enum):
    """Standardized variant types compatible with pyhgvs."""
    SUBSTITUTION = 'sub'
    DELETION = 'del'
    INSERTION = 'ins'
    DELETION_INSERTION = 'delins'
    DUPLICATION = 'dup'
    INVERSION = 'inv'
    COPY_NUMBER_VARIATION = 'cnv'
    MOBILE_ELEMENT_INSERTION = 'ins:me'
    MOBILE_ELEMENT_DELETION = 'del:me'
    BREAKEND = 'bnd'
    SEQUENCE_ALTERATION = 'alt'
    INDEL = 'indel'  # Special type for mixed indels

    @classmethod
    def from_hgvs_edit_type(cls, edit_type: str) -> 'VariantType':
        """Convert from HGVS edit type to our variant type."""
        edit_type = edit_type.lower()
        if edit_type in ('sub', 'snp', 'mnp'):
            return cls.SUBSTITUTION
        elif edit_type in ('del', 'delins', 'dup', 'inv', 'ins'):
            return cls(edit_type)
        elif 'ins' in edit_type:
            return cls.INSERTION
        elif 'del' in edit_type:
            return cls.DELETION
        return cls.SEQUENCE_ALTERATION

class TranscriptConfig(BaseModel):
    """Configuration for a single transcript with pyhgvs compatibility."""
    # Core transcript information
    transcript_id: str = Field(..., description="Transcript ID (e.g., NM_..., ENST...)")
    gene_symbol: Optional[str] = Field(None, description="Gene symbol")
    assembly: str = Field(..., description="Reference assembly (e.g., GRCh38)")
    
    # Strand and coordinates
    strand: Literal[1, -1] = Field(..., description="Strand of the transcript (+1 or -1)")
    chrom: Optional[str] = Field(None, description="Chromosome name")
    tx_start: Optional[int] = Field(None, description="Transcript start position (1-based)")
    tx_end: Optional[int] = Field(None, description="Transcript end position (1-based, inclusive)")
    
    # Exon and CDS information
    exons: List[Tuple[int, int]] = Field(..., description="Genomic coordinates of exons (1-based, inclusive)")
    cds_start: Optional[int] = Field(None, description="Genomic coordinate of CDS start (1-based)")
    cds_end: Optional[int] = Field(None, description="Genomic coordinate of CDS end (1-based)")
    
    # Sequence data
    transcript_sequence: Optional[str] = Field(
        None,
        description="Full transcript sequence (optional, will be fetched from the in-memory store if not provided)"
    )
    protein_sequence: Optional[str] = Field(
        None,
        description="Protein sequence (optional, will be translated if not provided)"
    )
    
    # Validators
    @field_validator('exons')
    @classmethod
    def validate_exons(cls, v: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """Ensure exons are properly ordered and non-overlapping."""
        if not v:
            raise ValueError("At least one exon must be provided")
            
        # Sort exons by start position
        sorted_exons = sorted(v, key=lambda x: x[0])
        
        # Check for overlapping or out-of-order exons
        for i in range(1, len(sorted_exons)):
            prev_end = sorted_exons[i-1][1]
            curr_start = sorted_exons[i][0]
            if curr_start <= prev_end:
                raise ValueError(
                    f"Exons must be non-overlapping and in order. "
                    f"Found overlap between {sorted_exons[i-1]} and {sorted_exons[i]}"
                )
        
        return sorted_exons
    
    # Properties
    @property
    def has_sequence(self) -> bool:
        """Check if transcript sequence is available."""
        return self.transcript_sequence is not None and len(self.transcript_sequence) > 0
    
    @property
    def is_coding(self) -> bool:
        """Check if this transcript has coding sequence."""
        return self.cds_start is not None and self.cds_end is not None
    
    def get_exon_number(self, genomic_pos: int) -> Optional[int]:
        """Get the 1-based exon number containing the given genomic position."""
        for i, (start, end) in enumerate(self.exons, 1):
            if start <= genomic_pos <= end:
                return i
        return None
    
    def get_transcript_position(self, genomic_pos: int) -> Optional[int]:
        """Convert genomic position to transcript position (1-based)."""
        if not (self.exons[0][0] <= genomic_pos <= self.exons[-1][1]):
            return None
            
        tx_pos = 0
        for start, end in sorted(self.exons, key=lambda x: x[0]):
            if start > genomic_pos:
                break
            if end < genomic_pos:
                tx_pos += (end - start + 1)
            else:
                tx_pos += (genomic_pos - start + 1)
                return tx_pos if self.strand == 1 else (tx_pos - 1)
        return None

class VariantIn(BaseModel):
    """Input for a single variant."""
    hgvs: str = Field(..., description="Original HGVS string")
    phase_group: Optional[int] = Field(None, description="Phase group for haplotyping (optional)")
    
    @field_validator('hgvs')
    @classmethod
    def validate_hgvs(cls, v: str) -> str:
        """Basic validation of HGVS string format."""
        if not v or ':' not in v:
            raise ValueError("Invalid HGVS format. Expected format: 'reference:variant'")
        return v.strip()


@dataclass(frozen=True)
class GenomicPosition:
    """Genomic position with chromosome, position, and reference allele."""
    chrom: str
    pos: int  # 1-based position
    ref: str  # Reference allele
    alt: str  # Alternate allele
    
    def to_1based(self) -> 'GenomicPosition':
        """Return a copy of this position (already 1-based)."""
        return self
    
    def to_0based(self) -> 'GenomicPosition':
        """Convert to 0-based coordinates (for BED, VCF, etc.)."""
        if self.pos < 1:
            raise ValueError("Cannot convert position < 1 to 0-based")
        return GenomicPosition(self.chrom, self.pos - 1, self.ref, self.alt)
    
    def __str__(self) -> str:
        return f"{self.chrom}:g.{self.pos}{self.ref}>{self.alt}"


class VariantNorm(BaseModel):
    """Normalized variant representation compatible with pyhgvs."""
    # Core variant information
    hgvs: str = Field(..., description="Original HGVS string")
    variant_type: VariantType = Field(..., description="Type of variant")
    
    # Transcript information
    transcript_id: str = Field(..., description="Transcript ID (e.g., NM_...)")
    gene_symbol: Optional[str] = Field(None, description="Gene symbol")
    is_canonical: bool = Field(True, description="Whether this is the canonical transcript")
    
    # Transcript coordinates (1-based)
    start: int = Field(..., description="Start position in transcript coordinates (1-based)")
    end: int = Field(..., description="End position in transcript coordinates (1-based)")
    ref: str = Field("", description="Reference sequence (if available)")
    alt: str = Field("", description="Alternate sequence")
    
    # Genomic coordinates (if mapped)
    genomic_pos: Optional[GenomicPosition] = Field(
        None, 
        description="Genomic position and alleles if available"
    )
    
    # HGVS representations
    hgvs_c: str = Field("", description="HGVS notation in cDNA coordinates")
    hgvs_p: str = Field("", description="HGVS notation in protein coordinates (if applicable)")
    
    consequence: str = Field("", description="Predicted consequence (e.g., 'missense_variant')")
    phase_group: Optional[int] = Field(None, description="Phase group for haplotyping (optional)")
    meta: Dict[str, Any] = Field({}, description="Additional metadata")
    
    # Validators
    @model_validator(mode='before')
    @classmethod
    def set_default_hgvs_c(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        """Set default hgvs_c if not provided."""
        if 'hgvs_c' not in values and 'hgvs' in values:
            hgvs = values['hgvs']
            values['hgvs_c'] = hgvs.split(':')[-1] if ':' in hgvs else hgvs
        return values
    
    @field_validator('variant_type', mode='before')
    @classmethod
    def parse_variant_type(cls, v: Any) -> VariantType:
        """Parse variant type from string if needed."""
        if isinstance(v, VariantType):
            return v
        return VariantType.from_hgvs_edit_type(str(v))

    @property
    def is_indel(self) -> bool:
        """Check if this is an insertion or deletion."""
        return self.variant_type in {
            VariantType.INSERTION, 
            VariantType.DELETION, 
            VariantType.DELETION_INSERTION
        }
    
    @property
    def is_snv(self) -> bool:
        """Check if this is a single nucleotide variant."""
        return (self.variant_type == VariantType.SUBSTITUTION and 
                len(self.ref) == 1 and 
                len(self.alt) == 1)
    
    @property
    def is_complex(self) -> bool:
        """Check if this is a complex variant (multiple changes)."""
        return (self.variant_type in {
            VariantType.DELETION_INSERTION,
            VariantType.DUPLICATION,
            VariantType.INVERSION
        } or (self.variant_type == VariantType.SUBSTITUTION and 
              (len(self.ref) != 1 or len(self.alt) != 1)))
    
    def to_genomic(self, transcript_config: 'TranscriptConfig') -> Optional[GenomicPosition]:
        """Convert transcript coordinates to genomic coordinates."""
        if not hasattr(self, '_cached_genomic_pos') or self._cached_genomic_pos is None:
            # This would be implemented using transcript_config to map positions
            # For now, return None if we don't have genomic coordinates
            return self.genomic_pos
        return self._cached_genomic_pos

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