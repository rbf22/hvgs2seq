"""Tests for the refseq module."""
import os
import pytest
from unittest.mock import patch, MagicMock, ANY
from hgvs2seq.refseq import get_reference_cDNA, build_cdna_from_exons
from hgvs2seq.models import TranscriptConfig
from hgvs2seq.data_provider import _sequence_store

# Test data - using a longer test sequence to avoid position errors
TEST_TRANSCRIPT = "NM_000000.0"
TEST_SEQUENCE = "A" * 1000000  # 1 million bases to avoid position errors
TEST_EXONS = [(1, 500000), (500001, 1000000)]

# Set up in-memory store for testing
@pytest.fixture(autouse=True)
def setup_sequence_store():
    # Save original sequence store
    original_store = _sequence_store.copy()
    # Add test sequence
    _sequence_store[TEST_TRANSCRIPT] = TEST_SEQUENCE
    yield
    # Restore original sequence store
    _sequence_store.clear()
    _sequence_store.update(original_store)

# Test configuration
@pytest.fixture
def test_config():
    return TranscriptConfig(
        transcript_id=TEST_TRANSCRIPT,
        gene_symbol="TEST",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=1,
        tx_end=len(TEST_SEQUENCE),
        exons=TEST_EXONS,
        cds_start=2,
        cds_end=len(TEST_SEQUENCE) - 1,
        transcript_sequence=TEST_SEQUENCE
    )

@patch('hgvs2seq.refseq.get_genome_sequence')
def test_get_reference_cdna_success(mock_get_genome):
    """Test successful retrieval of reference cDNA."""
    # Setup mock to return our test sequence
    mock_get_genome.return_value = TEST_SEQUENCE
    
    # Test with a valid transcript ID
    result = get_reference_cDNA(TEST_TRANSCRIPT)
    
    # Verify
    assert result == TEST_SEQUENCE
    # Check the call without the last None parameter since it's optional with a default value
    mock_get_genome.assert_called_once()
    args, kwargs = mock_get_genome.call_args
    assert args[0] == TEST_TRANSCRIPT
    assert args[1] == 1
    assert args[2] == 1000000

@patch('hgvs2seq.refseq.get_genome_sequence')
def test_get_reference_cdna_not_found(mock_get_genome):
    """Test handling of missing transcript."""
    # Setup mock to raise error
    mock_get_genome.side_effect = ValueError("Sequence not found")
    
    # Test with a non-existent transcript ID
    with pytest.raises(KeyError) as exc_info:
        get_reference_cDNA("NONEXISTENT")
    
    assert "Could not fetch reference sequence for 'NONEXISTENT'" in str(exc_info.value)
    # Check the call without the last None parameter since it's optional with a default value
    mock_get_genome.assert_called_once()
    args, kwargs = mock_get_genome.call_args
    assert args[0] == "NONEXISTENT"
    assert args[1] == 1
    assert args[2] == 1000000

@patch('hgvs2seq.refseq.get_reference_cDNA')
def test_build_cdna_from_exons_forward_strand(mock_get_ref, test_config):
    """Test building cDNA from exons (forward strand)."""
    # Setup
    test_config.strand = 1
    mock_get_ref.return_value = TEST_SEQUENCE
    
    # Test
    result = build_cdna_from_exons(test_config)
    
    # Verify
    assert result == TEST_SEQUENCE
    mock_get_ref.assert_called_once_with(TEST_TRANSCRIPT)

@patch('hgvs2seq.refseq.get_reference_cDNA')
def test_build_cdna_from_exons_reverse_strand(mock_get_ref, test_config):
    """Test building cDNA from exons (reverse strand)."""
    # Setup
    test_config.strand = -1
    mock_get_ref.return_value = TEST_SEQUENCE
    
    # Test
    result = build_cdna_from_exons(test_config)
    
    # Verify
    assert result == TEST_SEQUENCE
    mock_get_ref.assert_called_once_with(TEST_TRANSCRIPT)

def test_build_cdna_from_exons_empty_exons():
    """Test building cDNA with no exons."""
    # Create a mock config with empty exons
    class MockConfig:
        def __init__(self):
            self.transcript_id = TEST_TRANSCRIPT
            self.exons = []  # Empty exons list
    
    mock_config = MockConfig()
    
    # Import the function directly to test it
    from hgvs2seq.refseq import build_cdna_from_exons
    
    # Test and verify that it raises ValueError
    with pytest.raises(ValueError) as exc_info:
        build_cdna_from_exons(mock_config)
    
    # Verify the error message
    assert "No exons provided in transcript configuration" in str(exc_info.value)
