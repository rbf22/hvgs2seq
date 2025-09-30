"""Tests for the data_provider module."""
import sys
import os
import pytest
import logging
from unittest.mock import patch, MagicMock

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the data_provider module
import hgvs2seq.data_provider as dp

# Clear any existing cache
dp._transcript_cache = None

def test_get_transcript_data():
    """Test getting transcript data from the in-memory store."""
    # Reset the cache
    dp._transcript_cache = None
    
    # Test getting a transcript
    transcript = dp.get_transcript_data("NM_000000.0")
    
    # Verify the transcript data
    assert transcript is not None
    assert transcript['transcript_id'] == "NM_000000.0"
    assert transcript['gene_symbol'] == 'MOCK_GENE'
    assert transcript['chrom'] == 'chr1'
    assert transcript['strand'] == '+1'
    assert transcript['cds_start'] == 1000
    assert transcript['cds_end'] == 2000
    assert len(transcript['exons']) == 3
    assert len(transcript['transcript_sequence']) == 4000  # 1000 * 4 (ATCG)
    
    # Test that the transcript is cached
    with patch.dict(dp._transcript_cache, clear=True):
        dp._transcript_cache["NM_000000.0"] = {
            'transcript_id': "NM_000000.0",
            'gene_symbol': 'CACHED_GENE',
            'chrom': 'chr2',
            'strand': '-1',
            'cds_start': 2000,
            'cds_end': 3000,
            'exons': [(1, 200), (301, 400)],
            'transcript_sequence': 'ATCG' * 500,
        }
        cached_transcript = dp.get_transcript_data("NM_000000.0")
        assert cached_transcript['gene_symbol'] == 'CACHED_GENE'
        assert len(dp._transcript_cache) == 1


def test_get_genome_sequence():
    """Test getting a genomic sequence from the in-memory store."""
    # Test with a sequence that exists in the in-memory store
    result = dp.get_genome_sequence("chr1", 1, 10)
    assert result == 'N' * 10, "Unexpected sequence returned for chr1"
    
    # Test with a sequence that doesn't exist
    with pytest.raises(ValueError, match="Sequence ID not found: nonexistent_chr"):
        dp.get_genome_sequence("nonexistent_chr", 1, 10)
    
    # Test with invalid positions
    with pytest.raises(ValueError, match="Invalid positions: start=0, end=10 for sequence of length 10000"):
        dp.get_genome_sequence("chr1", 0, 10)  # start < 1
    
    with pytest.raises(ValueError, match="Invalid positions: start=1, end=10001 for sequence of length 10000"):
        dp.get_genome_sequence("chr1", 1, 10001)  # end > len(seq)
    
    with pytest.raises(ValueError, match="Invalid positions: start=10, end=5 for sequence of length 10000"):
        dp.get_genome_sequence("chr1", 10, 5)  # start > end


def test_get_genome_sequence_with_mock():
    """Test getting a genomic sequence with a mock sequence store."""
    # Create a test sequence
    test_seq = "ACGT" * 25  # 100 bases
    
    # Patch the _sequence_store for this test
    with patch.dict(dp._sequence_store, {'test_chr': test_seq}, clear=True):
        # Test getting a subsequence
        result = dp.get_genome_sequence("test_chr", 1, 4)
        assert result == "ACGT", "Unexpected sequence returned"
        
        # Test getting the entire sequence
        result = dp.get_genome_sequence("test_chr", 1, 100)
        assert result == test_seq, "Unexpected full sequence returned"
        
        # Test getting the last few bases
        result = dp.get_genome_sequence("test_chr", 97, 100)
        assert result == "ACGT"[-4:], "Unexpected end sequence returned"


def test_transcript_sequence_retrieval():
    """Test retrieving transcript sequences from the in-memory store."""
    # Test getting a transcript sequence
    result = dp.get_genome_sequence("NM_000000.0", 1, 4)
    assert result == "ATGC"[:4], "Unexpected transcript sequence returned"
    
    # Test getting a longer transcript sequence
    result = dp.get_genome_sequence("NM_000000.0", 1, 10)
    assert result == ("ATGC" * 3)[:10], "Unexpected transcript sequence returned"


def test_error_handling():
    """Test error handling in the data provider functions."""
    # Test with invalid sequence ID in get_genome_sequence
    with pytest.raises(ValueError, match="Sequence ID not found: invalid_id"):
        dp.get_genome_sequence("invalid_id", 1, 10)
    
    # Test with non-string seq_id - should raise ValueError from _sequence_store.get()
    with pytest.raises(ValueError, match="Sequence ID not found: 123"):
        dp.get_genome_sequence(123, 1, 10)  # seq_id is not str
    
    # Test with invalid position values in get_genome_sequence
    with pytest.raises(ValueError, match="Invalid positions: start=0, end=10 for sequence of length 10000"):
        dp.get_genome_sequence("chr1", 0, 10)  # start < 1
        
    with pytest.raises(ValueError, match="Invalid positions: start=1, end=10001 for sequence of length 10000"):
        dp.get_genome_sequence("chr1", 1, 10001)  # end > len(seq)
        
    with pytest.raises(ValueError, match="Invalid positions: start=10, end=5 for sequence of length 10000"):
        dp.get_genome_sequence("chr1", 10, 5)  # start > end
