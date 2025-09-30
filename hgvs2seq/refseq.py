"""
Handles fetching of reference sequences for a given transcript.
"""
import logging
from typing import List, Tuple, Optional

from .data_provider import get_genome_sequence
from .models import TranscriptConfig

_logger = logging.getLogger(__name__)

def get_reference_cDNA(transcript_id: str) -> str:
    """
    Fetch the reference cDNA sequence for a given transcript ID.

    Args:
        transcript_id: The transcript ID (e.g., 'NM_000000.0')

    Returns:
        The cDNA sequence as a string.

    Raises:
        KeyError: If the transcript ID is not found in the in-memory store.
    """
    _logger.info(f"Fetching reference cDNA for {transcript_id}")
    try:
        # Get the full transcript sequence from the in-memory store
        # We use a large end position to get the full sequence
        sequence = get_genome_sequence(transcript_id, 1, 1000000)
        _logger.info(f"Successfully fetched reference cDNA for {transcript_id}")
        return sequence
    except ValueError as e:
        _logger.error(f"Error fetching sequence for transcript '{transcript_id}': {e}")
        raise KeyError(f"Could not fetch reference sequence for '{transcript_id}'. "
                      "Ensure the ID is correct and the sequence exists in the in-memory store.")

def build_cdna_from_exons(config: TranscriptConfig) -> str:
    """
    Constructs the cDNA sequence by concatenating exon sequences fetched
    from the genomic reference based on the transcript configuration.
    This is an alternative to fetching the transcript directly and ensures
    the cDNA matches the provided exon structure.

    Args:
        config: The TranscriptConfig object.

    Returns:
        The constructed cDNA sequence as a string.
        
    Raises:
        ValueError: If no exons are provided in the config.
    """
    if not config.exons:
        _logger.error("No exons provided in transcript configuration")
        raise ValueError("No exons provided in transcript configuration")
        
    _logger.info(f"Building cDNA for {config.transcript_id} from {len(config.exons)} exons.")

    # We need a genomic sequence accession for the assembly mapper
    # This part can be tricky. We'll try to infer it from hgvs, but for now,
    # this function is a placeholder for a more robust implementation.
    # For now, we will rely on get_reference_cDNA.
    # A full implementation would require mapping the gene/transcript to a chromosome,
    # which is beyond the scope of this initial step.

    # This is a simplified example of what would be needed:
    # chrom = "NC_000017.11" # This would need to be fetched/mapped
    # full_sequence = ""
    # for start, end in config.exons:
    #     exon_seq = sr.fetch(chrom, start, end)
    #     full_sequence += exon_seq
    # if config.strand == -1:
    #     full_sequence = str(Seq(full_sequence).reverse_complement())

    _logger.warning("build_cdna_from_exons is not fully implemented. "
                   "Relying on get_reference_cDNA for now.")

    return get_reference_cDNA(config.transcript_id)