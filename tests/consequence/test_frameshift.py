"""Tests for the FrameshiftChecker."""
import pytest
from hgvs2seq.consequence.frameshift import FrameshiftChecker
from hgvs2seq.models import TranscriptConfig

@pytest.fixture
def minimal_config():
    """Return a minimal valid TranscriptConfig for testing."""
    return TranscriptConfig(
        transcript_id="test_transcript",
        assembly="GRCh38",
        strand=1,
        exons=[(1, 1000)],
        chrom="chr1",
        tx_start=1,
        tx_end=1000,
        cds_end=1000
    )

def test_frameshift_detection(minimal_config):
    """Test that a frameshift is correctly detected."""
    # Insertion of 1 base (frameshift)
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGGTACAGT",  # +1 insertion
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MVSX"  # Different length protein
    )
    assert result == "frameshift_variant"

def test_no_frameshift_when_inframe_indel(minimal_config):
    """Test that no frameshift is reported for inframe indels."""
    # Insertion of 3 bases (no frameshift)
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGGGTACGT",  # +3 insertion (inframe)
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MGV"
    )
    assert result is None

def test_frameshift_with_protein_validation(minimal_config):
    """Test frameshift detection using protein sequence validation."""
    # Same length DNA but different protein lengths
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGGTACGA",  # Same length but different protein
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MVSX"  # Different length protein
    )
    assert result == "frameshift_variant"

def test_no_frameshift_when_same_protein_length(minimal_config):
    """Test no frameshift when protein lengths are the same."""
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGGTACGA",  # Same length, different protein
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MVS"  # Same length protein
    )
    assert result is None

def test_frameshift_deletion(minimal_config):
    """Test detection of frameshift caused by a deletion."""
    # Deletion causing frameshift
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGTACGT",  # -1 deletion causing frameshift
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MX"
    )
    assert result == "frameshift_variant"

def test_empty_sequence_handling(minimal_config):
    """Test behavior with empty sequences."""
    # Empty reference sequence
    result = FrameshiftChecker.check(
        ref_cdna="",
        edited_cdna="ATG",
        config=minimal_config,
        ref_protein="",
        edited_protein="M"
    )
    assert result is None
    
    # Empty edited sequence (should be a frameshift)
    result = FrameshiftChecker.check(
        ref_cdna="ATG",
        edited_cdna="",
        config=minimal_config,
        ref_protein="M",
        edited_protein=""
    )
    assert result == "frameshift_variant"

def test_edge_cases(minimal_config):
    """Test various edge cases."""
    # Very short sequences
    result = FrameshiftChecker.check(
        ref_cdna="ATG",
        edited_cdna="AT",
        config=minimal_config,
        ref_protein="M",
        edited_protein=""
    )
    assert result == "frameshift_variant"
    
    # No change
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGGTACGT",
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MVR"
    )
    assert result is None
    assert result is None
    
    # Empty edited sequence
    result = FrameshiftChecker.check(
        ref_cdna="ATG",
        edited_cdna="",
        config=minimal_config,
        ref_protein="M",
        edited_protein=""
    )
    assert result is not None
    assert result.consequence == "frameshift_variant"

def test_frameshift_with_identical_proteins(minimal_config):
    """Test that no frameshift is reported when protein sequences are identical."""
    # Even if the cDNA lengths differ, if the proteins are the same, it's not a frameshift
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",  # MVR
        edited_cdna="ATGGGTACGT",  # +3 insertion, but same protein
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MVR"
    )
    assert result is None

def test_frameshift_with_missing_protein_sequences(minimal_config):
    """Test that None is returned when protein sequences are not provided."""
    # Missing ref_protein
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGGTACGTG",  # +1 insertion
        config=minimal_config,
        ref_protein="",
        edited_protein="MVRX"
    )
    assert result is None
    
    # Missing edited_protein
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGGTACGTG",  # +1 insertion
        config=minimal_config,
        ref_protein="MVR",
        edited_protein=""
    )
    assert result is None

def test_frameshift_with_short_sequences(minimal_config):
    """Test behavior with very short sequences."""
    # Very short sequences
    result = FrameshiftChecker.check(
        ref_cdna="A",
        edited_cdna="AT",  # +1 insertion
        config=minimal_config,
        ref_protein="",
        edited_protein=""
    )
    assert result is None  # No protein sequence provided
    
    # Short sequences with proteins
    result = FrameshiftChecker.check(
        ref_cdna="ATG",  # M
        edited_cdna="AT",  # Incomplete codon
        config=minimal_config,
        ref_protein="M",
        edited_protein="X"
    )
    assert result is not None
    assert result.consequence == "frameshift_variant"

def test_deletion_frameshift(minimal_config):
    """Test detection of frameshift caused by deletion."""
    # Deletion causing frameshift: ATG GTA CGT (MVR) -> ATG TAC GT (MX)
    result = FrameshiftChecker.check(
        ref_cdna="ATGGTACGT",
        edited_cdna="ATGTACGT",  # -1 deletion causing frameshift
        config=minimal_config,
        ref_protein="MVR",
        edited_protein="MX"
    )
    assert result is not None
    assert result.consequence == "frameshift_variant"
    # The frameshift occurs at position 2, changing Val to Ter (stop)
    assert result.hgvs_p == "p.Val2Terfs"
    assert result.protein_sequence == "MX"
    assert result.meta["frameshift_position"] == 2
    assert result.meta["frameshift_type"] == "deletion"
    assert result.meta["length_change"] == 1
