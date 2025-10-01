"""Tests for the StartLossChecker."""
import pytest
from hgvs2seq.consequence.start_loss import StartLossChecker

class TestStartLossChecker:
    """Test cases for StartLossChecker."""
    
    def test_start_loss_detection(self, create_mock_config):
        """Test basic start loss detection."""
        config = create_mock_config(cds_start=1)
        # ATG (Met) -> CTG (Leu)
        result = StartLossChecker.check(
            ref_cdna="ATGGTACGT",  # MVR
            edited_cdna="CTGGTACGT",  # LVR
            config=config,
            ref_protein="MVR",
            edited_protein="LVR"
        )
        assert result == "start_lost"
    
    def test_start_loss_with_coordinates(self, create_mock_config):
        """Test start loss with genomic coordinates."""
        config = create_mock_config(cds_start=100)
        # Variant at position 100-102 changes ATG->CTG
        result = StartLossChecker.check(
            ref_cdna="ATGGTACGT",  # MVR
            edited_cdna="CTGGTACGT",  # LVR
            config=config,
            ref_protein="MVR",
            edited_protein="LVR",
            variant_start=100,
            variant_end=102
        )
        assert result == "start_lost"
    
    def test_no_start_loss_when_no_change(self, create_mock_config):
        """Test no start loss when start codon is unchanged."""
        config = create_mock_config(cds_start=1)
        result = StartLossChecker.check(
            ref_cdna="ATGGTACGT",  # MVR
            edited_cdna="ATGGTACGA",  # MVR (silent mutation)
            config=config,
            ref_protein="MVR",
            edited_protein="MVR"
        )
        assert result is None
    
    def test_no_start_loss_when_not_atg(self, create_mock_config):
        """Test no start loss when reference doesn't start with ATG."""
        config = create_mock_config(cds_start=4)  # Start codon is at position 4
        
        result = StartLossChecker.check(
            ref_cdna="XXXATGGTACGT",  # XMVR
            edited_cdna="XXXCTGTACGT",  # XLR
            config=config,
            ref_protein="XMVR",
            edited_protein="XLR"
        )
        assert result is None
    
    def test_start_loss_to_stop(self, create_mock_config):
        """Test detection when start codon is changed to a stop codon."""
        config = create_mock_config(cds_start=1)
        # ATG (Met) -> TAG (Stop)
        result = StartLossChecker.check(
            ref_cdna="ATGGTACGT",  # MVR
            edited_cdna="TAGGTACGT",  # *VR
            config=config,
            ref_protein="MVR",
            edited_protein="*VR"
        )
        assert result == "start_lost"
    
    def test_short_sequence(self, create_mock_config):
        """Test behavior with sequences that are too short."""
        config = create_mock_config(cds_start=1)
        result = StartLossChecker.check(
            ref_cdna="AT",
            edited_cdna="CT",
            config=config,
            ref_protein="",
            edited_protein=""
        )
        assert result is None
    
    def test_no_cds_start_in_config(self):
        """Test behavior when config has no cds_start attribute."""
        # Create a config without cds_start
        class DummyConfig:
            pass
            
        config = DummyConfig()
        # Test case where start codon is lost (M->L)
        result = StartLossChecker.check(
            ref_cdna="ATGGTACGT",
            edited_cdna="CTGGTACGT",
            config=config,
            ref_protein="MVR",
            edited_protein="LVR"
        )
        # Should detect start loss because we can't verify coordinates
        assert result == "start_lost"
        
        # Test case where start codon is preserved
        result = StartLossChecker.check(
            ref_cdna="ATGGTACGT",
            edited_cdna="ATGGTACGA",
            config=config,
            ref_protein="MVR",
            edited_protein="MVR"
        )
        # Should return None since the start codon is preserved
        assert result is None
    
    def test_empty_sequences(self, create_mock_config):
        """Test behavior with empty sequences."""
        config = create_mock_config(cds_start=1)
        
        # Empty reference
        result = StartLossChecker.check(
            ref_cdna="",
            edited_cdna="ATGGTACGT",
            config=config,
            ref_protein="",
            edited_protein="MVR"
        )
        assert result is None
        
        # Empty edited
        result = StartLossChecker.check(
            ref_cdna="ATGGTACGT",
            edited_cdna="",
            config=config,
            ref_protein="MVR",
            edited_protein=""
        )
        assert result is None
