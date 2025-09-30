"""
Tests for FASTA output generation in hgvs2seq.io.fasta.
"""
import pytest
import sys
import os

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hgvs2seq.io.fasta import generate_fasta_output, _format_fasta_entry
from hgvs2seq.models import (
    TranscriptConfig, 
    TranscriptOutcome,
    ProteinOutcome,
    NMDOutcome,
    SequenceBundle,
    VariantType
)

# Test data
TEST_PROTEIN_SEQ = "MREK*"
TEST_CDNA_SEQ = "ATGCGTGAGAAATAG"

def test_format_fasta_entry():
    """Test the _format_fasta_entry helper function."""
    header = "test_header"
    sequence = "ATGCGT"
    expected = ">test_header\nATGCGT\n"
    assert _format_fasta_entry(header, sequence) == expected

def test_generate_fasta_output():
    """Test generating FASTA output from a SequenceBundle."""
    # Create a test config
    config = TranscriptConfig(
        transcript_id="NM_TEST.1",
        gene_symbol="TEST",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=1,
        tx_end=100,
        cds_start=10,
        cds_end=90,
        exons=[[1, 50], [60, 100]],
        transcript_sequence=TEST_CDNA_SEQ
    )
    
    # Create a test protein outcome
    protein_outcome = ProteinOutcome(
        protein_sequence=TEST_PROTEIN_SEQ,
        consequence="missense_variant"
    )
    
    # Create a test NMD outcome
    nmd_outcome = NMDOutcome(
        status="not_applicable",
        rationale="No premature termination codon"
    )
    
    # Create a test transcript outcome
    outcome = TranscriptOutcome(
        haplotype_id=1,
        scenario_id="test_scenario_1",
        mrna_sequence=TEST_CDNA_SEQ,
        protein_outcome=protein_outcome,
        nmd=nmd_outcome
    )
    
    # Create a sequence bundle
    bundle = SequenceBundle(
        primary_outcomes=[outcome],
        alternate_outcomes=[],
        provenance={
            'version': 'test_version',
            'timestamp': '2025-01-01T00:00:00Z',
            'source': 'test_source'
        }
    )
    
    # Generate FASTA output
    fasta_output = generate_fasta_output(bundle, config)
    
    # Check if the output contains expected sequences
    assert TEST_PROTEIN_SEQ in fasta_output
    assert TEST_CDNA_SEQ in fasta_output
    assert "type=protein" in fasta_output
    assert "type=cDNA" in fasta_output

def test_pah_sequence_loading():
    """Test that we can load and process the PAH sequence correctly."""
    # This test assumes the PAH sequence is available in the data directory
    from pathlib import Path
    
    # Path to the PAH sequence file (update this path as needed)
    pah_fasta = Path("data/GRCh38.chr12.fa")
    
    if not pah_fasta.exists():
        pytest.skip(f"PAH sequence file not found at {pah_fasta}")
    
    # Create a simple FASTA loader function for testing
    def load_fasta(fasta_path: str) -> str:
        """Simple FASTA loader for testing."""
        with open(fasta_path, 'r') as f:
            # Skip header
            next(f)
            # Read sequence
            return ''.join(line.strip() for line in f if not line.startswith('>'))
    
    # Load the sequence
    sequence = load_fasta(str(pah_fasta))
    
    # Basic checks
    assert len(sequence) > 0, "Sequence should not be empty"
    
    # Check sequence composition (should be mostly A, C, G, T, N)
    bases = set(sequence.upper())
    valid_bases = {'A', 'C', 'G', 'T', 'N'}
    assert bases.issubset(valid_bases), f"Unexpected bases in sequence: {bases - valid_bases}"
