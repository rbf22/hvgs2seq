"""
Tests for the NMD (nonsense-mediated decay) prediction module.
"""
import pytest
from hgvs2seq.models import ProteinOutcome, TranscriptConfig, NMDOutcome
from hgvs2seq.consequence.nmd import check_nmd, NMD_RULE_DISTANCE_NT

def test_check_nmd_not_applicable_no_ptc():
    """Test NMD check when there's no premature termination codon."""
    # Create a protein outcome without a PTC
    protein_outcome = ProteinOutcome(
        consequence="missense",
        protein_sequence="MADEV*"
    )
    
    # Create a dummy transcript config (shouldn't matter for this case)
    config = TranscriptConfig(
        transcript_id="NM_123456.7",
        gene_symbol="TEST1",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=100,
        tx_end=400,
        exons=[(100, 200), (300, 400)],
        cds_start=150,
        cds_end=350
    )
    
    # Check NMD
    result = check_nmd(protein_outcome, config)
    
    # Verify the result
    assert result.status == "not_applicable"
    assert "No premature termination codon" in result.rationale

def test_check_nmd_intronless_transcript():
    """Test NMD check with an intronless transcript."""
    # Create a protein outcome with a PTC
    protein_outcome = ProteinOutcome(
        consequence="stop_gain",
        protein_sequence="MA*DEV"  # PTC at position 2 (0-based)
    )
    
    # Create a transcript config with a single exon (intronless)
    config = TranscriptConfig(
        transcript_id="NM_123456.7",
        gene_symbol="TEST1",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=100,
        tx_end=200,
        exons=[(100, 200)],  # Single exon
        cds_start=110,
        cds_end=170
    )
    
    # Check NMD
    result = check_nmd(protein_outcome, config)
    
    # Verify the result
    assert result.status == "escape"
    assert "Transcript is intronless" in result.rationale

def test_check_nmd_ptc_in_last_exon():
    """Test NMD check when PTC is in the last exon."""
    # Create a protein outcome with a PTC
    protein_outcome = ProteinOutcome(
        consequence="stop_gain",
        protein_sequence="MADE*V"  # PTC at position 4 (0-based)
    )
    
    # Create a transcript config with two exons
    # PTC will be at c.122 (cds_start + 4*3 = 110 + 12 = 122)
    # Last junction is at 100 (length of first exon)
    # So PTC is after the last junction (in last exon)
    config = TranscriptConfig(
        transcript_id="NM_123456.7",
        gene_symbol="TEST1",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=100,
        tx_end=300,
        exons=[(100, 199), (200, 300)],  # First exon: 100nt, Second exon: 101nt
        cds_start=110,  # CDS starts at 110 in the first exon
        cds_end=250     # CDS ends at 250 in the second exon
    )
    
    # Check NMD
    result = check_nmd(protein_outcome, config)
    
    # Verify the result
    assert result.status == "escape"
    assert "PTC at c.122 is located in the last exon" in result.rationale

def test_check_nmd_ptc_before_last_junction_far():
    """Test NMD check when PTC is far from the last exon-exon junction."""
    # Create a protein outcome with a PTC
    protein_outcome = ProteinOutcome(
        consequence="stop_gain",
        protein_sequence="M*ADEV"  # PTC at position 1 (0-based)
    )
    
    # Create a transcript config with two exons
    # PTC will be at c.13 (cds_start + 1*3 = 10 + 3 = 13)
    # Last junction is at 100 (length of first exon)
    # Distance to last junction is 87nt (100 - 13)
    # NMD_RULE_DISTANCE_NT is 55, so this should trigger NMD
    config = TranscriptConfig(
        transcript_id="NM_123456.7",
        gene_symbol="TEST1",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=10,
        tx_end=300,
        exons=[(100, 199), (200, 300)],  # First exon: 100nt, Second exon: 101nt
        cds_start=10,  # CDS starts at 10 in the first exon
        cds_end=250    # CDS ends at 250 in the second exon
    )
    
    # Check NMD
    result = check_nmd(protein_outcome, config)
    
    # Verify the result
    assert result.status == "likely"
    assert f"PTC at c.13 is {100-13} nt upstream of the last exon-exon junction" in result.rationale
    assert f"which is >= {NMD_RULE_DISTANCE_NT} nt" in result.rationale

def test_check_nmd_ptc_before_last_junction_near():
    """Test NMD check when PTC is near the last exon-exon junction."""
    # Create a protein outcome with a PTC
    protein_outcome = ProteinOutcome(
        consequence="stop_gain",
        protein_sequence="MADE*V"  # PTC at position 4 (0-based)
    )
    
    # Create a transcript config with two exons
    # PTC will be at c.13 (cds_start + 4*3 = 1 + 12 = 13)
    # Last junction is at 20 (length of first exon)
    # Distance to last junction is 7nt (20 - 13)
    # NMD_RULE_DISTANCE_NT is 55, so this should not trigger NMD
    config = TranscriptConfig(
        transcript_id="NM_123456.7",
        gene_symbol="TEST1",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=1,
        tx_end=300,
        exons=[(100, 119), (200, 300)],  # First exon: 20nt, Second exon: 101nt
        cds_start=1,  # CDS starts at 1 in the first exon
        cds_end=250   # CDS ends at 250 in the second exon
    )
    
    # Check NMD
    result = check_nmd(protein_outcome, config)
    
    # Verify the result
    assert result.status == "escape"
    assert f"PTC at c.13 is {20-13} nt upstream of the last exon-exon junction" in result.rationale
    assert f"which is < {NMD_RULE_DISTANCE_NT} nt" in result.rationale

def test_check_nmd_no_cds_start():
    """Test NMD check when CDS start is not defined."""
    # Create a protein outcome with a PTC
    protein_outcome = ProteinOutcome(
        consequence="stop_gain",
        protein_sequence="MA*DEV"  # PTC at position 2 (0-based)
    )
    
    # Create a transcript config without CDS start but with multiple exons
    config = TranscriptConfig(
        transcript_id="NM_123456.7",
        gene_symbol="TEST1",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=100,
        tx_end=300,
        exons=[(100, 150), (200, 300)],  # Two exons to avoid the intronless check
        cds_start=None,  # No CDS start defined
        cds_end=280
    )
    
    # Check NMD
    result = check_nmd(protein_outcome, config)
    
    # Verify the result
    assert result.status == "not_applicable"
    assert "CDS start is not defined" in result.rationale
