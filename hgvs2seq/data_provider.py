"""
Simple in-memory data provider for testing purposes.
This provides a minimal implementation of sequence and transcript data access
without requiring external dependencies.
"""
import logging
import os
from typing import Dict, Any, Optional, Type, TypeVar, Union, TypeAlias

_logger = logging.getLogger(__name__)

# Flag to track if SeqRepo is available
SEQREPO_AVAILABLE = False

# In-memory sequence store for testing
_sequence_store = {
    'chr1': 'N' * 10000,  # Dummy chromosome sequence
    'NM_000000.0': 'ATGC' * 1000,  # Dummy transcript sequence
}

# Global variables to hold the initialized connections (singletons)
_transcript_cache = None

# Type aliases
TranscriptData = Dict[str, Any]

# Type variable for sequence store
try:
    from typing_extensions import TypeAlias
except ImportError:
    # For Python < 3.10
    from typing import Type as TypeAlias  # type: ignore

# Dummy SeqRepo-like class for type hints
class DummySeqRepo:
    """Dummy SeqRepo implementation for type hints and testing."""
    def get_sequence(self, seq_id: str, start: int, end: int) -> str:
        """Get a sequence slice."""
        seq = _sequence_store.get(seq_id, '')
        return seq[start-1:end]  # Convert to 0-based for slicing

# Type aliases for compatibility
SeqRepo: TypeAlias = DummySeqRepo
SeqRepoType: TypeAlias = DummySeqRepo


def get_transcript_data(transcript_id: str) -> Optional[TranscriptData]:
    """
    Get transcript data for the given transcript ID.
    
    This is a simplified implementation that would typically query a database
    or use a local cache of transcript information.
    
    Args:
        transcript_id: The transcript ID (e.g., 'NM_000000.0')
        
    Returns:
        Dict containing transcript information or None if not found
    """
    global _transcript_cache
    
    # Initialize the cache if it doesn't exist
    if _transcript_cache is None:
        _transcript_cache = {}
        # In a real implementation, you would load this from a database or file
        # For now, we'll just return a mock transcript
        _logger.warning("Using mock transcript data. In production, load from a real data source.")
    
    # Check if we have this transcript in the cache
    if transcript_id not in _transcript_cache:
        # In a real implementation, you would fetch this from a database
        # For now, we'll return a mock transcript
        _transcript_cache[transcript_id] = {
            'transcript_id': transcript_id,
            'gene_symbol': 'MOCK_GENE',
            'chrom': 'chr1',
            'strand': '+1',
            'cds_start': 1000,
            'cds_end': 2000,
            'exons': [(1, 100), (201, 300), (401, 500)],
            'transcript_sequence': 'ATCG' * 1000,  # Mock sequence
        }
    
    return _transcript_cache.get(transcript_id)


def get_sr(seqrepo_class=None):
    """
    Returns a singleton SeqRepo instance.
    
    Args:
        seqrepo_class: Optional class to use for creating the SeqRepo instance.
                      If not provided, uses DummySeqRepo.
                      
    Returns:
        A SeqRepo instance.
        
    Raises:
        EnvironmentError: If SEQREPO_ROOT is not set or the path does not exist.
    """
    global _sr
    
    if _sr is not None:
        return _sr
        
    # Check if SEQREPO_ROOT is set
    seqrepo_root = os.environ.get('SEQREPO_ROOT')
    if not seqrepo_root:
        raise EnvironmentError("SEQREPO_ROOT environment variable is not set")
        
    # Check if the path exists
    if not os.path.exists(seqrepo_root):
        raise EnvironmentError(f"SEQREPO_ROOT path does not exist: {seqrepo_root}")
    
    # Use the provided class or default to DummySeqRepo
    sr_class = seqrepo_class or DummySeqRepo
    _sr = sr_class(seqrepo_root)
    _logger.info(f"Initialized {sr_class.__name__} with root: {seqrepo_root}")
    
    return _sr


def get_genome_sequence(seq_id: str, start: int, end: int, sr_instance=None):
    """
    Get a genomic sequence.
    
    Args:
        seq_id: The sequence ID (e.g., 'chr1' or 'NM_000000.0')
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)
        sr_instance: Optional SeqRepo instance to use. If not provided, uses the in-memory store.
        
    Returns:
        The genomic sequence as a string
        
    Raises:
        ValueError: If the sequence is not found or positions are invalid
    """
    try:
        if sr_instance is not None:
            # Use the provided SeqRepo instance
            # Note: SeqRepo uses 0-based, half-open intervals [start, end)
            return sr_instance.fetch(seq_id, start - 1, end)
        else:
            # Fall back to in-memory store
            seq = _sequence_store.get(seq_id)
            if seq is None:
                raise ValueError(f"Sequence ID not found: {seq_id}")
                
            if start < 1 or end > len(seq) or start > end:
                raise ValueError(
                    f"Invalid positions: start={start}, end={end} for sequence of length {len(seq)}"
                )
                
            return seq[start-1:end]  # Convert to 0-based for slicing
    except Exception as e:
        error_msg = f"Error getting sequence {seq_id}:{start}-{end}: {str(e)}"
        _logger.error(error_msg)
        raise ValueError(error_msg) from e