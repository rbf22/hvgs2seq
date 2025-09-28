import pytest
from unittest.mock import MagicMock
from hgvs2seq.parse import parse_and_normalize_variants
from hgvs2seq.models import VariantIn, TranscriptConfig

# Define more accurate mock objects that reflect the hgvs library structure
class MockHgvsVariant:
    def __init__(self, ac, type, posedit):
        self.ac = ac
        self.type = type
        self.posedit = posedit
    def __str__(self):
        return f"{self.ac}:{self.type}.{self.posedit}"

class MockPosEdit:
    def __init__(self, pos, edit):
        self.pos = pos
        self.edit = edit
    def __str__(self):
        return f"{self.pos}{self.edit}"

class MockInterval:
    """Mocks the 'pos' attribute, which has 'start' and 'end'."""
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __str__(self):
        return f"{self.start}_{self.end}"

class MockPosition:
    """Mocks the 'start' and 'end' position objects."""
    def __init__(self, base, offset=0):
        self.base = base
        self.offset = offset
    def __str__(self):
        if self.offset == 0:
            return str(self.base)
        return f"{self.base}{'+' if self.offset > 0 else ''}{self.offset}"

class MockEdit:
    def __init__(self, type, alt=""):
        self.type = type
        self.alt = alt
    def __str__(self):
        return f"delins{self.alt}" if self.type == 'sub' else f"{self.type}{self.alt}"


@pytest.fixture
def mock_data_providers(monkeypatch):
    """Mocks the hgvs data providers (AssemblyMapper and UTA) and the hgvs parser."""
    mock_am = MagicMock()
    mock_hdp = MagicMock()
    mock_parser = MagicMock()

    # The normalize_variant method should just return the variant it was given in the mock
    mock_hdp.normalize_variant.side_effect = lambda v: v

    # Define the behavior of the mocked parser
    def parse_side_effect(hgvs_string):
        if "c.15G>T" in hgvs_string:
            pos = MockInterval(MockPosition(15), MockPosition(15))
            edit = MockEdit("sub", "T")
            return MockHgvsVariant("NM_TOY001.1", "c", MockPosEdit(pos, edit))
        if "c.50+2T>G" in hgvs_string:
            pos = MockInterval(MockPosition(50, 2), MockPosition(50, 2))
            edit = MockEdit("sub", "G")
            return MockHgvsVariant("NM_TOY001.1", "c", MockPosEdit(pos, edit))
        return MagicMock()

    mock_parser.parse_hgvs_variant.side_effect = parse_side_effect

    # Patch the getter functions and the parser instance
    monkeypatch.setattr("hgvs2seq.parse.get_am", lambda: mock_am)
    monkeypatch.setattr("hgvs2seq.parse.get_hdp", lambda: mock_hdp)
    monkeypatch.setattr("hgvs2seq.parse._parser", mock_parser)

@pytest.fixture
def transcript_config():
    from hgvs2seq.config import load_config
    return load_config("tests/fixtures/test_transcript_config.json")

def test_parse_simple_variant(mock_data_providers, transcript_config):
    """Tests parsing a simple cDNA variant."""
    variants_in = [VariantIn(hgvs="NM_TOY001.1:c.15G>T", phase_group=1)]
    result = parse_and_normalize_variants(variants_in, transcript_config)

    assert len(result) == 1
    norm_var = result[0]
    assert norm_var.kind == "sub"
    assert norm_var.c_start == 15
    assert norm_var.c_end == 15
    assert norm_var.alt == "T"
    assert norm_var.phase_group == 1
    assert norm_var.c_start_offset == 0

def test_parse_splice_site_variant(mock_data_providers, transcript_config):
    """Tests parsing a variant with an intronic offset."""
    variants_in = [VariantIn(hgvs="NM_TOY001.1:c.50+2T>G")]
    result = parse_and_normalize_variants(variants_in, transcript_config)

    assert len(result) == 1
    norm_var = result[0]
    assert norm_var.kind == "sub"
    assert norm_var.c_start == 50
    assert norm_var.c_start_offset == 2
    assert norm_var.alt == "G"