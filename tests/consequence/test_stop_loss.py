"""Tests for the StopLossChecker."""
import pytest
from hgvs2seq.consequence.stop_loss import StopLossChecker
from tests.conftest import create_mock_config

class TestStopLossChecker:
    """Test cases for StopLossChecker."""
    
    def test_stop_loss_basic(self, create_mock_config):
        """Test basic stop loss detection."""
        config = create_mock_config()
        # TAA (stop) -> CAA (Q)
        result = StopLossChecker.check(
            ref_cdna="ATGTAA",  # M*
            edited_cdna="ATGCAA",  # MQ
            config=config,
            ref_protein="M*",
            edited_protein="MQ",
            variant_start=4,
            variant_end=6
        )
        assert result == "stop_lost"
    
    def test_stop_loss_with_coordinates(self, create_mock_config):
        """Test stop loss with genomic coordinates."""
        config = create_mock_config(cds_start=100, cds_end=105)  # ATGTAA (M*) at 100-105
        # Variant at position 103-105 changes TAA->CAA (*->Q)
        result = StopLossChecker.check(
            ref_cdna="ATGTAA",  # M*
            edited_cdna="ATGCAA",  # MQ
            config=config,
            ref_protein="M*",
            edited_protein="MQ",
            variant_start=103,
            variant_end=105
        )
        # The checker should detect stop loss regardless of coordinates for this test
        assert result == "stop_lost"
    
    def test_no_stop_loss_wrong_position(self, create_mock_config):
        """Test no stop loss when variant is not at the stop codon."""
        config = create_mock_config(cds_start=100, cds_end=105)  # ATGTAA (M*) at 100-105
        # Variant at position 1-3 changes ATG->CTG (M->L)
        result = StopLossChecker.check(
            ref_cdna="ATGTAA",  # M*
            edited_cdna="CTGTAA",  # L*
            config=config,
            ref_protein="M*",
            edited_protein="L*",
            variant_start=1,
            variant_end=3
        )
        assert result is None
    
    def test_no_stop_loss_no_stop_in_reference(self, create_mock_config):
        """Test no stop loss when there's no stop in the reference."""
        config = create_mock_config()
        # No stop in reference or edited
        result = StopLossChecker.check(
            ref_cdna="ATGCAG",  # MQ
            edited_cdna="ATGTAG",  # M*
            config=config,
            ref_protein="MQ",
            edited_protein="M*"
        )
        assert result is None
