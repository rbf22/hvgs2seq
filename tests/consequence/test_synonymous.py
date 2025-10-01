"""Tests for the SynonymousChecker."""
import pytest
from hgvs2seq.consequence.synonymous import SynonymousChecker
from tests.conftest import create_mock_config

class TestSynonymousChecker:
    """Test cases for SynonymousChecker."""
    
    def test_synonymous_basic(self, create_mock_config):
        """Test basic synonymous variant detection."""
        config = create_mock_config()
        # CCT (P) -> CCC (P) - both code for Proline
        result = SynonymousChecker.check(
            ref_cdna="ATGCCT",  # MP
            edited_cdna="ATGCCC",  # MP
            config=config,
            ref_protein="MP",
            edited_protein="MP"
        )
        assert result == "synonymous_variant"
    
    def test_not_synonymous_different_protein(self, create_mock_config):
        """Test that non-synonymous variants are not marked as synonymous."""
        config = create_mock_config()
        # CCT (P) -> ACT (T) - different amino acid
        result = SynonymousChecker.check(
            ref_cdna="ATGCCT",  # MP
            edited_cdna="ATGACT",  # MT
            config=config,
            ref_protein="MP",
            edited_protein="MT"
        )
        assert result is None
    
    def test_not_synonymous_same_sequence(self, create_mock_config):
        """Test that identical sequences are not marked as synonymous."""
        config = create_mock_config()
        # No change
        result = SynonymousChecker.check(
            ref_cdna="ATGCCT",  # MP
            edited_cdna="ATGCCT",  # MP
            config=config,
            ref_protein="MP",
            edited_protein="MP"
        )
        assert result is None
    
    def test_not_synonymous_empty_sequences(self, create_mock_config):
        """Test behavior with empty sequences."""
        config = create_mock_config()
        # Empty sequences
        result = SynonymousChecker.check(
            ref_cdna="",
            edited_cdna="",
            config=config,
            ref_protein="",
            edited_protein=""
        )
        assert result is None
