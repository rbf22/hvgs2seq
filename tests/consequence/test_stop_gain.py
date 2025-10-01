"""Tests for the StopGainChecker."""
import pytest
from hgvs2seq.consequence.stop_gain import StopGainChecker
from tests.conftest import create_mock_config

class TestStopGainChecker:
    """Test cases for StopGainChecker."""
    
    def test_stop_gain_basic(self, create_mock_config):
        """Test basic stop gain detection."""
        config = create_mock_config()
        # CAG (Q) -> TAG (stop)
        result = StopGainChecker.check(
            ref_cdna="ATGCAG",  # MQ
            edited_cdna="ATGTAG",  # M*
            config=config,
            ref_protein="MQ",
            edited_protein="M*",
            variant_start=4,
            variant_end=6
        )
        assert result == "stop_gained"
    
    def test_stop_gain_with_coordinates(self, create_mock_config):
        """Test stop gain with genomic coordinates."""
        config = create_mock_config(cds_start=100)
        # Variant at position 100-102 changes CAG->TAG (Q->*)
        result = StopGainChecker.check(
            ref_cdna="ATGCAG",  # MQ
            edited_cdna="ATGTAG",  # M*
            config=config,
            ref_protein="MQ",
            edited_protein="M*",
            variant_start=100,
            variant_end=102
        )
        # The checker should detect stop gain regardless of coordinates for this test
        assert result == "stop_gained"
    
    def test_no_stop_gain_wrong_position(self, create_mock_config):
        """Test no stop gain when variant is not at the right position."""
        config = create_mock_config(cds_start=100)
        # Variant at position 1-3 changes ATG->CTG (M->L)
        result = StopGainChecker.check(
            ref_cdna="ATGCAG",  # MQ
            edited_cdna="CTGCAG",  # LQ
            config=config,
            ref_protein="MQ",
            edited_protein="LQ",
            variant_start=1,
            variant_end=3
        )
        assert result is None
    
    def test_no_stop_gain_natural_stop(self, create_mock_config):
        """Test no stop gain when the stop is already in the reference."""
        config = create_mock_config()
        # TAG (stop) -> TAA (stop) - both are stop codons
        result = StopGainChecker.check(
            ref_cdna="ATGTAG",  # M*
            edited_cdna="ATGTAA",  # M*
            config=config,
            ref_protein="M*",
            edited_protein="M*"
        )
        assert result is None
