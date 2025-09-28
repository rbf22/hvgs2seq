import json
from typing import TextIO
from ..models import SequenceBundle

def write_json(bundle: SequenceBundle, file_handle: TextIO, indent: int = 2):
    """
    Serializes a SequenceBundle to a JSON file.

    Args:
        bundle: The SequenceBundle to serialize.
        file_handle: A writable file handle.
        indent: The indentation level for the JSON output.
    """
    # pydantic's model_dump() is the modern way to get a serializable dict
    bundle_dict = bundle.model_dump(mode="json")
    json.dump(bundle_dict, file_handle, indent=indent)