import pytest
import json
from hgvs2seq.config import load_config
from hgvs2seq.models import TranscriptConfig

# Path to the valid test fixture
VALID_CONFIG_PATH = "tests/fixtures/test_transcript_config.json"

def test_load_valid_config():
    """
    Tests that a valid transcript configuration file is loaded correctly.
    """
    # Load the configuration
    config = load_config(VALID_CONFIG_PATH)

    # Assert that the returned object is of the correct type
    assert isinstance(config, TranscriptConfig)

    # Load the raw data to compare values
    with open(VALID_CONFIG_PATH, 'r') as f:
        raw_data = json.load(f)

    # Assert that the attributes match the file content
    assert config.transcript_id == raw_data["transcript_id"]
    assert config.gene_symbol == raw_data["gene_symbol"]
    assert config.assembly == raw_data["assembly"]
    assert config.strand == raw_data["strand"]
    assert config.exons == [tuple(e) for e in raw_data["exons"]]
    assert config.cds_start_c == raw_data["cds_start_c"]
    assert config.cds_end_c == raw_data["cds_end_c"]

def test_load_nonexistent_config():
    """
    Tests that loading a non-existent file raises FileNotFoundError.
    """
    with pytest.raises(FileNotFoundError):
        load_config("tests/fixtures/nonexistent_file.json")

def test_load_invalid_json(tmp_path):
    """
    Tests that loading a file with invalid JSON raises a ValueError.
    """
    invalid_json_file = tmp_path / "invalid.json"
    invalid_json_file.write_text("{'transcript_id': 'NM_123', 'assembly': 'GRCh38'") # Missing closing brace

    with pytest.raises(ValueError, match="Invalid JSON"):
        load_config(str(invalid_json_file))

def test_load_config_with_validation_error(tmp_path):
    """
    Tests that loading a config with missing required fields raises a ValueError.
    """
    invalid_data_file = tmp_path / "invalid_data.json"
    # Missing 'transcript_id' which is a required field
    invalid_data_file.write_text('{"assembly": "GRCh38", "strand": 1, "exons": [[1, 100]]}')

    with pytest.raises(ValueError, match="Configuration validation failed"):
        load_config(str(invalid_data_file))