"""
Handles generation of JSON formatted outputs.
"""
import logging
from ..models import SequenceBundle

_logger = logging.getLogger(__name__)

def generate_json_output(bundle: SequenceBundle, indent: int = 2) -> str:
    """
    Generates a JSON formatted string from the SequenceBundle.

    Args:
        bundle: The SequenceBundle object containing all results.
        indent: The indentation level for pretty-printing the JSON.

    Returns:
        A JSON formatted string.
    """
    _logger.info("Generating JSON output from SequenceBundle.")
    try:
        # Pydantic models have a built-in method to serialize to a JSON string.
        # `model_dump_json` is the method for Pydantic v2.
        return bundle.model_dump_json(indent=indent)
    except Exception as e:
        _logger.error(f"Failed to serialize SequenceBundle to JSON: {e}")
        raise TypeError(f"Error during JSON serialization: {e}") from e