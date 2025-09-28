import json
import yaml
import pytest
from hgvs2seq.config import TranscriptConfig, load_config

@pytest.fixture
def json_config_file(tmp_path):
    """Creates a temporary JSON config file for testing."""
    config_data = {
        "transcript_id": "NM_000059.3",
        "gene_symbol": "BRCA1",
        "assembly": "GRCh38",
        "strand": -1,
        "exons": [(1, 100), (201, 300)],
        "cds_start_c": 50,
        "cds_end_c": 250,
    }
    file_path = tmp_path / "config.json"
    with open(file_path, "w") as f:
        json.dump(config_data, f)
    return str(file_path)

@pytest.fixture
def yaml_config_file(tmp_path):
    """Creates a temporary YAML config file for testing."""
    config_data = {
        "transcript_id": "NM_000059.3",
        "gene_symbol": "BRCA1",
        "assembly": "GRCh38",
        "strand": -1,
        "exons": [[1, 100], [201, 300]],
        "cds_start_c": 50,
        "cds_end_c": 250,
    }
    file_path = tmp_path / "config.yaml"
    with open(file_path, "w") as f:
        yaml.dump(config_data, f)
    return str(file_path)

def test_load_config_from_json(json_config_file):
    """Tests loading a TranscriptConfig from a JSON file."""
    config = load_config(json_config_file)
    assert isinstance(config, TranscriptConfig)
    assert config.transcript_id == "NM_000059.3"
    assert config.gene_symbol == "BRCA1"
    assert config.assembly == "GRCh38"
    assert config.strand == -1
    assert config.exons == [(1, 100), (201, 300)]
    assert config.cds_start_c == 50
    assert config.cds_end_c == 250

def test_load_config_from_yaml(yaml_config_file):
    """Tests loading a TranscriptConfig from a YAML file."""
    config = load_config(yaml_config_file)
    assert isinstance(config, TranscriptConfig)
    assert config.transcript_id == "NM_000059.3"
    assert config.exons == [(1, 100), (201, 300)] # Note: YAML loads tuples as lists

def test_load_config_invalid_file_type(tmp_path):
    """Tests that a ValueError is raised for unsupported file types."""
    # Create a dummy file with an unsupported extension
    invalid_file = tmp_path / "config.txt"
    invalid_file.write_text("content")

    with pytest.raises(ValueError, match="must be a .json or .yaml file"):
        load_config(str(invalid_file))