import pytest
from hgvs2seq.consequence import ConsequenceAnalyzer
from hgvs2seq.models import TranscriptConfig

@pytest.fixture
def analyzer():
    return ConsequenceAnalyzer()

@pytest.fixture
def minimal_config():
    """Return a minimal valid TranscriptConfig for testing."""
    return TranscriptConfig(
        transcript_id="test_transcript",
        assembly="GRCh38",
        strand=1,
        exons=[(1, 1000)],  # (start, end) tuples for exons
        chrom="chr1",
        tx_start=1,
        tx_end=1000,
        cds_start=1,
        cds_end=1000
    )

def test_start_loss(analyzer, minimal_config):
    """Test detection of start loss variant."""
    ref_cdna = "ATGGTACGT"  # MVR
    edited_cdna = "CTGGTACGT"  # LV* (if it were translated)
    
    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)
    assert result.consequence == "start_lost"

def test_stop_gain(analyzer, minimal_config):
    """Test detection of stop gain variant."""
    ref_cdna = "ATGGTACGT"  # MVR
    # Change to TAA (stop) at position 4-6 (0-based 3-5)
    edited_cdna = "ATGTAA" + "CGT"  # M*R (premature stop)
    
    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)

def test_missense(analyzer, minimal_config):
    """Test detection of missense variant."""
    ref_cdna = "ATGGTACGT"  # MVR
    edited_cdna = "ATGGTCCGT"  # MVP

    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)
    # Check that the protein sequence is correctly translated
    assert result.protein_sequence == "MVP"
    # Check that missense_variant is in the list of consequences
    assert "missense_variant" in result.all_consequences

def test_synonymous(analyzer, minimal_config):
    """Test detection of synonymous variant."""
    ref_cdna = "ATGGTACGT"  # MVR
    edited_cdna = "ATGGTACGC"  # MVR (CGT and CGC both code for R)
    
    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)
    assert result.consequence == "synonymous_variant"

def test_frameshift(analyzer, minimal_config):
    """Test detection of frameshift variant."""
    ref_cdna = "ATGGTACGT"  # MVR
    edited_cdna = "ATGGTACGGT"  # MVR (with +1 insertion)

    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)
    # Check that the protein sequence is correctly translated
    # The first two codons should be correct, then frameshift
    assert result.protein_sequence.startswith("MVR")
    # Check that frameshift_variant is in the list of consequences
    assert "frameshift_variant" in result.all_consequences

def test_inframe_indel(analyzer, minimal_config):
    """Test detection of inframe insertion/deletion."""
    # Test inframe insertion
    ref_cdna = "ATGGTACGT"  # MVR
    edited_cdna = "ATGGTAAACGT"  # MVR (with +3 insertion)

    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)
    # Check that the protein sequence is correctly translated
    # The insertion should be in-frame, so the protein should be longer
    assert len(result.protein_sequence) > 3  # Original was 3 AAs
    # Check that inframe_insertion is in the list of consequences
    assert "inframe_insertion" in result.all_consequences
    
    # Test inframe deletion
    edited_cdna = "ATGGT"  # MV (with -3 deletion)
    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)
    assert result.consequence == "inframe_deletion"
def test_stop_loss(analyzer, minimal_config):
    """Test detection of stop loss variant."""
    # Create a sequence with a stop codon
    ref_cdna = "ATGGTATAA"  # MV*
    # Change stop to a different amino acid
    edited_cdna = "ATGGTACAA"  # MVQ
    
    result = analyzer.analyze(ref_cdna, edited_cdna, minimal_config)
    assert result.consequence == "stop_lost"
