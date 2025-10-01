"""Tests for the ConsequenceAnalyzer."""
import pytest
from hgvs2seq.consequence.analyzer import ConsequenceAnalyzer
from tests.conftest import create_mock_config

class TestAnalyzer:
    """Test cases for ConsequenceAnalyzer."""
    
    def test_analyze_stop_gain(self, create_mock_config):
        """Test analyzer with a stop gain variant."""
        analyzer = ConsequenceAnalyzer()
        config = create_mock_config(cds_start=100)
        
        # CAG (Q) -> TAG (stop)
        result = analyzer.analyze(
            ref_cdna="ATGCAG",  # MQ
            edited_cdna="ATGTAG",  # M*
            config=config,
            ref_protein="MQ",
            edited_protein="M*"
        )
        assert result == "stop_gained"
    
    def test_analyze_stop_loss(self, create_mock_config):
        """Test analyzer with a stop loss variant."""
        analyzer = ConsequenceAnalyzer()
        config = create_mock_config(cds_start=100)
        
        # TAA (stop) -> CAA (Q)
        result = analyzer.analyze(
            ref_cdna="ATGTAA",  # M*
            edited_cdna="ATGCAA",  # MQ
            config=config,
            ref_protein="M*",
            edited_protein="MQ"
        )
        assert result == "stop_lost"
    
    def test_analyze_synonymous(self, create_mock_config):
        """Test analyzer with a synonymous variant."""
        analyzer = ConsequenceAnalyzer()
        config = create_mock_config()
        
        # CCT (P) -> CCC (P) - both code for Proline
        result = analyzer.analyze(
            ref_cdna="ATGCCT",  # MP
            edited_cdna="ATGCCC",  # MP
            config=config,
            ref_protein="MP",
            edited_protein="MP"
        )
        assert result == "synonymous_variant"
    
    def test_analyze_missense(self, create_mock_config):
        """Test analyzer with a missense variant."""
        analyzer = ConsequenceAnalyzer()
        config = create_mock_config()
        
        # CCT (P) -> ACT (T)
        result = analyzer.analyze(
            ref_cdna="ATGCCT",  # MP
            edited_cdna="ATGACT",  # MT
            config=config,
            ref_protein="MP",
            edited_protein="MT"
        )
        assert result == "missense_variant"
    
    def test_analyze_inframe_indel(self, create_mock_config):
        """Test analyzer with an inframe insertion."""
        analyzer = ConsequenceAnalyzer()
        config = create_mock_config()
        
        # Insertion of 3 bases (inframe)
        # Original: ATG GTA CGT (MVR)
        # Insertion: ATG GGG TAC GT (MGVC)
        result = analyzer.analyze(
            ref_cdna="ATGGTACGT",  # MVR
            edited_cdna="ATGGGGTACGT",  # MGVC
            config=config,
            ref_protein="MVR",
            edited_protein="MGVC"
        )
        # This is actually a frameshift because the length difference is 2 (not a multiple of 3)
        assert result == "frameshift_variant"
        
        # Test a real inframe insertion (length difference is a multiple of 3)
        result = analyzer.analyze(
            ref_cdna="ATGGTACGT",  # MVR
            edited_cdna="ATGGGGGTACGT",  # MGVRC
            config=config,
            ref_protein="MVR",
            edited_protein="MGVR"
        )
        assert result == "inframe_insertion"
    
    def test_analyze_frameshift(self, create_mock_config):
        """Test analyzer with a frameshift variant."""
        analyzer = ConsequenceAnalyzer()
        config = create_mock_config()
        
        # Insertion of 1 base (frameshift)
        result = analyzer.analyze(
            ref_cdna="ATGGTACGT",  # MVR
            edited_cdna="ATGGGTACGT",  # MVSX
            config=config,
            ref_protein="MVR",
            edited_protein="MVSX"
        )
        assert result == "frameshift_variant"
