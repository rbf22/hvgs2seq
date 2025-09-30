"""
Tests for the data models in hgvs2seq.models.
"""
import pytest
from hgvs2seq.models import (
    TranscriptConfig, 
    VariantNorm, 
    VariantType,
    ProteinOutcome,
    NMDOutcome,
    TranscriptOutcome,
    SequenceBundle,
    GenomicPosition
)

# Test data for TranscriptConfig
transcript_config_data = {
    "transcript_id": "NM_123456.7",
    "gene_symbol": "TEST1",
    "assembly": "GRCh38",
    "strand": 1,
    "chrom": "chr1",
    "tx_start": 1000,
    "tx_end": 2000,
    "cds_start": 1100,
    "cds_end": 1900,
    "exons": [(1000, 1200), (1500, 1700), (1800, 2000)],
    "transcript_sequence": "ATGCGTACGT" * 100,
    "protein_sequence": "MREK*"
}

# Test data for GenomicPosition
genomic_pos_data = {
    "chrom": "chr1",
    "pos": 123456,
    "ref": "A",
    "alt": "T"
}

# Test data for VariantNorm
variant_norm_data = {
    "hgvs": "NM_123456.7:c.100A>T",
    "variant_type": VariantType.SUBSTITUTION,
    "transcript_id": "NM_123456.7",
    "start": 100,
    "end": 100,
    "ref": "A",
    "alt": "T",
    "hgvs_c": "NM_123456.7:c.100A>T",
    "hgvs_p": "p.Lys34Ter",
    "consequence": "stop_gain",
    "phase_group": None,
    "genomic_pos": GenomicPosition(**genomic_pos_data)
}

# Test data for ProteinOutcome
protein_outcome_data = {
    "protein_sequence": "MREK*",
    "hgvs_p": "p.Lys34Ter",
    "consequence": "stop_gain"
}

# Test data for NMDOutcome
nmd_outcome_data = {
    "status": "likely",
    "rationale": "Premature stop codon in last exon"
}

# Test data for TranscriptOutcome
transcript_outcome_data = {
    "haplotype_id": 0,
    "scenario_id": "baseline",
    "mrna_sequence": "ATGCGTACGT" * 100,
    "protein_outcome": ProteinOutcome(**protein_outcome_data),
    "nmd": NMDOutcome(**nmd_outcome_data),
    "confidence": 0.95
}

# Test data for SequenceBundle
sequence_bundle_data = {
    "primary_outcomes": [TranscriptOutcome(**transcript_outcome_data)],
    "alternate_outcomes": [],
    "provenance": {"version": "1.0.0", "tool": "hgvs2seq"},
    "warnings": ["No warnings"]
}

def test_transcript_config_creation():
    """Test creation of TranscriptConfig with valid data."""
    config = TranscriptConfig(**transcript_config_data)
    assert config.transcript_id == "NM_123456.7"
    assert config.gene_symbol == "TEST1"
    assert config.assembly == "GRCh38"
    assert config.strand == 1
    assert config.chrom == "chr1"
    assert config.tx_start == 1000
    assert config.tx_end == 2000
    assert config.cds_start == 1100
    assert config.cds_end == 1900
    assert config.exons == [(1000, 1200), (1500, 1700), (1800, 2000)]
    assert len(config.transcript_sequence) == 1000
    assert config.protein_sequence == "MREK*"

def test_transcript_config_optional_fields():
    """Test creation of TranscriptConfig with optional fields."""
    data = transcript_config_data.copy()
    del data["gene_symbol"]
    del data["chrom"]
    del data["tx_start"]
    del data["tx_end"]
    
    config = TranscriptConfig(**data)
    assert config.gene_symbol is None
    assert config.chrom is None
    assert config.tx_start is None
    assert config.tx_end is None

def test_variant_norm_creation():
    """Test creation of VariantNorm with valid data."""
    variant = VariantNorm(**variant_norm_data)
    assert variant.hgvs == "NM_123456.7:c.100A>T"
    assert variant.transcript_id == "NM_123456.7"
    assert variant.start == 100
    assert variant.end == 100
    assert variant.ref == "A"
    assert variant.alt == "T"
    assert variant.hgvs_c == "NM_123456.7:c.100A>T"
    assert variant.hgvs_p == "p.Lys34Ter"
    assert variant.consequence == "stop_gain"
    assert variant.phase_group is None

def test_variant_norm_optional_fields():
    """Test creation of VariantNorm with optional fields."""
    data = variant_norm_data.copy()
    del data["hgvs_p"]
    del data["phase_group"]
    
    variant = VariantNorm(**data)
    assert variant.hgvs_p == ""  # Default is empty string, not None
    assert variant.phase_group is None

def test_protein_outcome_creation():
    """Test creation of ProteinOutcome."""
    outcome = ProteinOutcome(protein_sequence="MREK*", hgvs_p="p.Lys34Ter", consequence="stop_gain")
    assert outcome.protein_sequence == "MREK*"
    assert outcome.consequence == "stop_gain"

def test_genomic_position():
    """Test GenomicPosition functionality."""
    pos = GenomicPosition(**genomic_pos_data)
    
    # Test 0-based conversion
    pos_0based = pos.to_0based()
    assert pos_0based.pos == 123455  # Should be 0-based now
    
    # Test 1-based conversion (should be idempotent)
    pos_1based = pos.to_1based()
    assert pos_1based.pos == 123456  # Should still be 1-based

def test_variant_norm_methods():
    """Test methods of VariantNorm class."""
    # Create a variant for testing
    variant = VariantNorm(**variant_norm_data)
    
    # Test properties
    assert not variant.is_indel  # This is a substitution
    assert variant.is_snv  # Single nucleotide variant
    assert not variant.is_complex  # Simple substitution

def test_protein_outcome():
    """Test ProteinOutcome model."""
    outcome = ProteinOutcome(**protein_outcome_data)
    assert outcome.protein_sequence == "MREK*"
    assert outcome.hgvs_p == "p.Lys34Ter"
    assert outcome.consequence == "stop_gain"

def test_nmd_outcome():
    """Test NMDOutcome model."""
    outcome = NMDOutcome(**nmd_outcome_data)
    assert outcome.status == "likely"
    assert "Premature stop" in outcome.rationale

def test_transcript_outcome():
    """Test TranscriptOutcome model."""
    outcome = TranscriptOutcome(**transcript_outcome_data)
    assert outcome.haplotype_id == 0
    assert outcome.scenario_id == "baseline"
    assert len(outcome.mrna_sequence) == 1000
    assert isinstance(outcome.protein_outcome, ProteinOutcome)
    assert isinstance(outcome.nmd, NMDOutcome)
    assert 0 <= outcome.confidence <= 1.0

def test_sequence_bundle():
    """Test SequenceBundle model."""
    bundle = SequenceBundle(**sequence_bundle_data)
    assert len(bundle.primary_outcomes) == 1
    assert len(bundle.alternate_outcomes) == 0
    assert "version" in bundle.provenance
    assert len(bundle.warnings) == 1
