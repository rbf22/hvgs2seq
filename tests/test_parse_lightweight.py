"""Lightweight tests for hgvs2seq.parse module without hgvs dependency."""
import pytest
from unittest.mock import MagicMock, patch
from hgvs2seq.models import VariantIn, TranscriptConfig, VariantType
from hgvs2seq.parse import (
    parse_and_normalize_variants,
    _project_genomic_to_transcript,
    _project_noncoding_to_coding,
    HGVSParseError,
    HGVSDataNotAvailableError
)

# Lightweight mock classes
class MockHgvsVariant:
    def __init__(self, ac, type, posedit, **kwargs):
        self.ac = ac
        self.type = type
        self.posedit = posedit
        self._data = kwargs
        
    def __str__(self):
        return f"{self.ac}:{self.type}.{self.posedit}"
        
    def __getattr__(self, name):
        if name in self._data:
            return self._data[name]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

class MockPosEdit:
    def __init__(self, pos, edit):
        self.pos = pos
        self.edit = edit
        
    def __str__(self):
        return f"{self.pos}{self.edit}"

class MockInterval:
    def __init__(self, start, end, uncertain=False):
        self.start = start
        self.end = end
        self.uncertain = uncertain
        
    def __str__(self):
        return f"{self.start}_{self.end}"

class MockPosition:
    def __init__(self, base, offset=0):
        self.base = base
        self.offset = offset
        
    def __str__(self):
        if self.offset == 0:
            return str(self.base)
        return f"{self.base}{'+' if self.offset > 0 else ''}{self.offset}"

class MockEdit:
    def __init__(self, type, ref=None, alt=""):
        self.type = type
        self.alt = alt
        self.ref = ref or ("" if type != "sub" else "A")
        
    def __str__(self):
        if self.type == 'sub':
            return f"{self.ref}>{self.alt}"
        elif self.type == 'del':
            return f"del{self.ref}"
        elif self.type == 'ins':
            return f"ins{self.alt}"
        elif self.type == 'delins':
            return f"del{self.ref}ins{self.alt}"
        elif self.type == 'dup':
            return f"dup{self.ref}"
        return f"{self.type}{self.alt}"

@pytest.fixture
def mock_parse_hgvs(monkeypatch):
    """Mock the parser.parse_hgvs_variant method."""
    def parse_hgvs(hgvs_string):
        if hgvs_string == "NM_000059.4:c.68A>G":
            pos = MockInterval(MockPosition(68), MockPosition(68))
            edit = MockEdit("sub", "A", "G")
            return MockHgvsVariant("NM_000059.4", "c", MockPosEdit(pos, edit))
        elif hgvs_string == "invalid_hgvs":
            raise HGVSParseError("Invalid HGVS format")
        return None
    
    mock = MagicMock()
    mock.parse_hgvs_variant.side_effect = parse_hgvs
    monkeypatch.setattr("hgvs2seq.parse._parser", mock)
    return mock

@pytest.fixture
def mock_transcript_data(monkeypatch):
    """Mock transcript data provider."""
    def get_transcript_data(transcript_id):
        return {
            'transcript_id': transcript_id,
            'exons': [(1, 1000)],
            'cds_start': 50,
            'cds_end': 950,
            'strand': 1,
            'transcript_sequence': 'A' * 1000
        }
    
    monkeypatch.setattr("hgvs2seq.parse.get_transcript_data", get_transcript_data)
    return get_transcript_data

def test_parse_simple_variant(mock_parse_hgvs, mock_transcript_data):
    """Test parsing a simple substitution variant."""
    variants_in = [VariantIn(hgvs="NM_000059.4:c.68A>G")]
    transcript_config = TranscriptConfig(
        transcript_id="NM_000059.4",
        gene_symbol="BRCA2",
        assembly="GRCh38",
        chrom="13",
        strand=1,
        tx_start=1,
        tx_end=1000,
        cds_start=50,
        cds_end=950,
        exons=[(1, 1000)]
    )
    
    result = parse_and_normalize_variants(variants_in, transcript_config)
    
    assert len(result) == 1
    variant = result[0]
    assert variant.transcript_id == "NM_000059.4"
    assert variant.variant_type == VariantType.SUBSTITUTION
    assert variant.hgvs == "NM_000059.4:c.68A>G"

def test_parse_invalid_hgvs(mock_parse_hgvs, mock_transcript_data):
    """Test handling of invalid HGVS strings."""
    # Test invalid HGVS format in VariantIn validation
    with pytest.raises(ValueError):
        VariantIn(hgvs="invalid_hgvs")
    
    # Test invalid HGVS that passes VariantIn but fails parsing
    variants_in = [VariantIn(hgvs="invalid_hgvs:variant")]
    transcript_config = TranscriptConfig(
        transcript_id="NM_000059.4",
        gene_symbol="BRCA2",
        assembly="GRCh38",
        chrom="13",
        strand=1,
        tx_start=1,
        tx_end=1000,
        cds_start=50,
        cds_end=950,
        exons=[(1, 1000)]
    )
    
    with pytest.raises(ValueError):
        parse_and_normalize_variants(variants_in, transcript_config)

def test_project_genomic_to_transcript():
    """Test projection of genomic to transcript coordinates."""
    mock_variant = MockHgvsVariant(
        "chr13", "g",
        MockPosEdit(
            MockInterval(MockPosition(32315640), MockPosition(32315640)),
            MockEdit("sub", "A", "G")
        )
    )
    
    transcript_data = {
        'transcript_id': 'NM_000059.4',
        'exons': [(32315474, 32315670), (32316520, 32316686)],
        'strand': 1,
        'cds_start': 32315640,
        'cds_end': 33344547
    }
    
    result = _project_genomic_to_transcript(mock_variant, transcript_data)
    assert result.ac == 'NM_000059.4'
    assert result.type == 'c'

def test_project_noncoding_to_coding():
    """Test projection of non-coding to coding coordinates."""
    mock_variant = MockHgvsVariant(
        "NR_123456.1", "n",
        MockPosEdit(
            MockInterval(MockPosition(100), MockPosition(100)),
            MockEdit("sub", "A", "G")
        )
    )
    
    transcript_data = {
        'transcript_id': 'NM_000059.4',
        'cds_start': 50,
        'cds_end': 200,
        'strand': 1,
        'exons': [(1, 200)],
        'transcript_sequence': 'A' * 200,
        'chrom': '13',
        'tx_start': 1,
        'tx_end': 200
    }
    
    # Test the function
    result = _project_noncoding_to_coding(mock_variant, transcript_data)
    assert result.ac == 'NR_123456.1'  # Should keep the original accession
    assert result.type == 'c'  # Should change type to 'c'
