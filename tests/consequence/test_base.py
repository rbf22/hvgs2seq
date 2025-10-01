"""Tests for the ConsequenceChecker base class."""
import pytest
from hgvs2seq.consequence.base import ConsequenceChecker

class TestTranslateSequence:
    """Tests for the translate_sequence method."""
    
    def test_translate_sequence_with_start_codon(self):
        """Test translation starting from a start codon."""
        # Test a simple sequence with ATG start codon
        dna = "XXXATGTTTAAACCC"  # ATG (M), TTT (F), AAA (K), CCC (P)
        expected = "MFP"  # Stops at first stop codon (TAA)
        assert ConsequenceChecker.translate_sequence(dna) == expected
        
    def test_translate_sequence_no_start_codon(self):
        """Test translation when no start codon is present."""
        # Sequence with no ATG
        dna = "TTTTTTAAACCC"
        assert ConsequenceChecker.translate_sequence(dna) == ""
        
    def test_translate_sequence_stop_codon(self):
        """Test that translation stops at the first stop codon."""
        # Sequence with stop codon in the middle
        dna = "ATGTGTTTATAAAGGG"  # ATG (M), TGT (C), TTA (L), TAA (stop), AGG (R)
        expected = "MCL"  # Should stop at TAA
        assert ConsequenceChecker.translate_sequence(dna) == expected
        
    def test_translate_sequence_short_sequence(self):
        """Test with sequences that are too short."""
        assert ConsequenceChecker.translate_sequence("") == ""
        assert ConsequenceChecker.translate_sequence("A") == ""
        assert ConsequenceChecker.translate_sequence("AT") == ""
        
    def test_translate_sequence_uppercase_lowercase(self):
        """Test that the method handles mixed case input."""
        dna_upper = "ATGTGT"  # M, C
        dna_lower = "atgtgt"  # M, C
        dna_mixed = "AtGtGt"  # M, C
        expected = "MC"
        
        assert ConsequenceChecker.translate_sequence(dna_upper) == expected
        assert ConsequenceChecker.translate_sequence(dna_lower) == expected
        assert ConsequenceChecker.translate_sequence(dna_mixed) == expected
        
    def test_translate_sequence_with_invalid_codons(self):
        """Test that invalid codons are translated to 'X'."""
        dna = "ATGNNNTAG"  # ATG (M), NNN (X), TAG (stop)
        assert ConsequenceChecker.translate_sequence(dna) == "MX"
