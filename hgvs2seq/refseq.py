import hgvs.dataproviders.uta
from functools import lru_cache

# This is a simplified setup. In a real application, this might be
# part of a shared data provider connection management system.
try:
    hdp = hgvs.dataproviders.uta.connect()
except Exception as e:
    raise RuntimeError(
        "Failed to connect to HGVS data provider (UTA) for sequence fetching. "
        "Please check network connection and UTA service status. "
        f"Original error: {e}"
    )

@lru_cache(maxsize=128)
def get_reference_sequence(transcript_id: str) -> str:
    """
    Fetches the reference cDNA sequence for a given transcript ID from UTA.

    Args:
        transcript_id: The RefSeq or Ensembl transcript ID (e.g., "NM_000551.3").

    Returns:
        The cDNA sequence as a string.

    Raises:
        ValueError: If the sequence could not be fetched.
    """
    try:
        # UTA provides sequence data via the `get_seq` method.
        # For transcripts, the accession is the ID itself.
        sequence = hdp.get_seq(transcript_id)
        if not sequence:
            raise ValueError(f"No sequence found for transcript '{transcript_id}'.")
        return sequence
    except Exception as e:
        # Catch other potential exceptions from the data provider.
        raise ValueError(
            f"An error occurred while fetching the sequence for '{transcript_id}': {e}"
        )