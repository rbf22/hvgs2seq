"""Tests for the project module."""
import pytest
from unittest.mock import MagicMock, patch, PropertyMock

import pytest

from hgvs2seq.hgvs_compat import SequenceVariant, Interval, SimplePosition, PosEdit
from hgvs2seq.project import project_variant, project_c_to_g
from hgvs2seq.models import TranscriptConfig, GenomicPosition, VariantNorm, VariantType

# Test data
TEST_TRANSCRIPT = "NM_000000.0"
TEST_GENE = "TEST"
TEST_ASSEMBLY = "GRCh38"
TEST_STRAND = 1
TEST_EXONS = [(1000, 2000), (3000, 4000)]
TEST_CDS_START = 1500
TEST_CDS_END = 3500
TEST_CHROM = "chr1"

@pytest.fixture
def test_config():
    """Create a test transcript configuration."""
    return TranscriptConfig(
        transcript_id=TEST_TRANSCRIPT,
        gene_symbol=TEST_GENE,
        assembly=TEST_ASSEMBLY,
        strand=TEST_STRAND,
        chrom=TEST_CHROM,
        tx_start=1000,
        tx_end=4000,
        exons=TEST_EXONS,
        cds_start=TEST_CDS_START,
        cds_end=TEST_CDS_END,
    )

@pytest.fixture
def mock_variant():
    """Create a mock VariantNorm object for testing."""
    variant = MagicMock(spec=VariantNorm)
    variant.hgvs = 'NM_000000.0:c.123A>G'
    variant.hgvs_c = 'c.123A>G'
    variant.transcript_id = TEST_TRANSCRIPT
    variant.variant_type = VariantType.SUBSTITUTION
    variant.start = 123
    variant.end = 123
    variant.ref = 'A'
    variant.alt = 'G'
    variant.genomic_pos = None
    
    # Mock the model_copy method
    variant.model_copy.return_value = variant
    
    return variant

def test_project_variant_same_transcript(mock_variant, test_config):
    """Test projecting a variant that's already in the target transcript."""
    # When the variant is already in the target transcript
    result = project_variant(mock_variant, test_config)
    
    # Should return the same variant
    assert result is mock_variant

def test_project_variant_c_to_n(mock_variant, test_config):
    """Test projecting a variant from c. to n. coordinates."""
    # Setup mock variant as cDNA variant
    mock_variant.hgvs = 'NM_000000.0:c.123A>G'
    mock_variant.hgvs_c = 'c.123A>G'
    
    # Call the function directly since we're not mocking anything
    result = project_variant(mock_variant, test_config)
    
    # Should return a variant (in this case, it should be the same variant)
    assert result is not None

def test_project_variant_unsupported_type(mock_variant, test_config):
    """Test that projecting a variant with an unsupported type is handled correctly."""
    # Setup mock variant with unsupported type
    mock_variant.hgvs = 'NM_000000.0:unsupported.123A>G'
    mock_variant.hgvs_c = 'unsupported.123A>G'
    
    # The current implementation doesn't raise an error for unsupported types,
    # it just returns the variant as-is. This test verifies that behavior.
    result = project_variant(mock_variant, test_config)
    assert result is mock_variant  # Should return the same variant object

def test_project_c_to_g_success(mock_variant, test_config):
    """Test successful projection from c. to g. coordinates."""
    # Setup mock variant with required attributes
    mock_variant.hgvs = 'NM_000000.0:c.123A>G'
    mock_variant.hgvs_c = 'c.123A>G'
    
    # Mock the model_copy method to return a copy with genomic position
    def mock_copy(update=None):
        if update:
            mock_variant.hgvs = update.get('hgvs', mock_variant.hgvs)
            mock_variant.genomic_pos = update.get('genomic_pos')
            mock_variant.variant_type = update.get('variant_type', mock_variant.variant_type)
        return mock_variant
    
    mock_variant.model_copy.side_effect = mock_copy
    
    # Test with config parameter
    result = project_c_to_g(mock_variant, config=test_config)
    assert result is not None
    assert result.genomic_pos is not None  # Should have genomic position

def test_project_c_to_g_failure():
    """Test handling of projection failure with invalid input."""
    # Test with invalid variant (None)
    with pytest.raises(ValueError, match="Variant cannot be None"):
        project_c_to_g(None)
    
    # Test with variant missing required attributes
    variant = MagicMock(spec=VariantNorm)
    variant.hgvs = 'NM_000000.0:invalid.123A>G'
    variant.hgvs_c = 'invalid.123A>G'
    variant.transcript_id = 'NM_000000.0'
    variant.genomic_pos = None
    with pytest.raises(ValueError, match="must be in c. coordinates"):
        project_c_to_g(variant)
