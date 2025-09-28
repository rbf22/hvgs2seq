import pytest
from hgvs2seq.parse import VariantNorm
from hgvs2seq.apply.single import apply_single_variant
from hgvs2seq.apply.batch import apply_batch

# A sample reference sequence for testing: ATG GCC TCA TGA TAG CAT CAT CAT CAT C
REF_SEQ = "ATGGCCTCATGATAGCATCATCATCATC"

@pytest.fixture
def ref_seq() -> str:
    """Provides the reference sequence for tests."""
    return REF_SEQ

# =============================================================================
# Tests for apply_single_variant
# =============================================================================

def test_apply_sub(ref_seq: str):
    """Test a single nucleotide substitution."""
    # c.4G>T changes the sequence from ATGGCCT... to ATGTCCT...
    variant = VariantNorm(hgvs_c="c.4G>T", kind="sub", c_start=4, c_end=4, alt="T", meta={})
    expected_seq = "ATGTCCTCATGATAGCATCATCATCATC"
    result_seq = apply_single_variant(ref_seq, variant)
    assert result_seq == expected_seq

def test_apply_del(ref_seq: str):
    """Test a deletion of a few bases."""
    # c.4_6del removes GCC from ATG-GCC-TCA... to become ATG-TCA...
    variant = VariantNorm(hgvs_c="c.4_6del", kind="del", c_start=4, c_end=6, alt=None, meta={})
    expected_seq = "ATGTCATGATAGCATCATCATCATC"
    result_seq = apply_single_variant(ref_seq, variant)
    assert result_seq == expected_seq

def test_apply_ins(ref_seq: str):
    """Test a simple insertion."""
    # c.6_7insA inserts 'A' between base 6 (C) and 7 (T)
    # ATG-GCC-TCA... becomes ATG-GCC-A-TCA...
    variant = VariantNorm(hgvs_c="c.6_7insA", kind="ins", c_start=6, c_end=7, alt="A", meta={})
    expected_seq = "ATGGCCATCATGATAGCATCATCATCATC"
    result_seq = apply_single_variant(ref_seq, variant)
    assert result_seq == expected_seq

def test_apply_delins(ref_seq: str):
    """Test a deletion followed by an insertion."""
    # c.4_6delinsTT removes GCC and inserts TT
    # ATG-GCC-TCA... becomes ATG-TT-TCA...
    variant = VariantNorm(hgvs_c="c.4_6delinsTT", kind="delins", c_start=4, c_end=6, alt="TT", meta={})
    expected_seq = "ATGTTTCATGATAGCATCATCATCATC"
    result_seq = apply_single_variant(ref_seq, variant)
    assert result_seq == expected_seq

def test_apply_mnv(ref_seq: str):
    """Test a multi-nucleotide variant (a delins of the same length)."""
    # c.4_6delinsTTA removes GCC and inserts TTA
    # ATG-GCC-TCA... becomes ATG-TTA-TCA...
    variant = VariantNorm(hgvs_c="c.4_6delinsTTA", kind="mnv", c_start=4, c_end=6, alt="TTA", meta={})
    expected_seq = "ATGTTATCATGATAGCATCATCATCATC"
    result_seq = apply_single_variant(ref_seq, variant)
    assert result_seq == expected_seq

def test_apply_dup(ref_seq: str):
    """Test a duplication."""
    # c.4_6dup duplicates GCC
    # ATG-GCC-TCA... becomes ATG-GCC-GCC-TCA...
    variant = VariantNorm(hgvs_c="c.4_6dup", kind="dup", c_start=4, c_end=6, alt=None, meta={})
    expected_seq = "ATGGCCGCCTCATGATAGCATCATCATCATC"
    result_seq = apply_single_variant(ref_seq, variant)
    assert result_seq == expected_seq

# =============================================================================
# Tests for apply_batch
# =============================================================================

def test_apply_batch_non_overlapping(ref_seq: str):
    """Test applying a batch of non-overlapping variants."""
    variants = [
        VariantNorm(hgvs_c="c.4_6del", kind="del", c_start=4, c_end=6, alt=None, meta={}),      # Deletes GCC
        VariantNorm(hgvs_c="c.10_11insA", kind="ins", c_start=10, c_end=11, alt="A", meta={}) # Inserts A between T and G
    ]
    # Starting sequence: ATGGCCTCATGATAGCATCATCATCATC
    # 1. Apply c.10_11insA: ATGGCCTCAT-A-GATAGCAT...
    # 2. Apply c.4_6del on the new sequence: ATG-TCAT-A-GATAGCAT...
    expected_seq = "ATGTCATAGATAGCATCATCATCATC"
    result_seq = apply_batch(ref_seq, variants)
    assert result_seq == expected_seq

def test_apply_batch_order_independence(ref_seq: str):
    """Test that the order of variants in the input list doesn't matter."""
    variants_shuffled = [
        VariantNorm(hgvs_c="c.10_11insA", kind="ins", c_start=10, c_end=11, alt="A", meta={}),
        VariantNorm(hgvs_c="c.4_6del", kind="del", c_start=4, c_end=6, alt=None, meta={})
    ]
    expected_seq = "ATGTCATAGATAGCATCATCATCATC"
    result_seq = apply_batch(ref_seq, variants_shuffled)
    assert result_seq == expected_seq

def test_out_of_bounds_error(ref_seq: str):
    """Test that a variant outside the sequence bounds raises an error."""
    variant = VariantNorm(hgvs_c="c.100G>A", kind="sub", c_start=100, c_end=100, alt="A", meta={})
    with pytest.raises(ValueError, match="out of bounds"):
        apply_single_variant(ref_seq, variant)

def test_unsupported_kind_error(ref_seq: str):
    """Test that an unsupported variant kind raises an error."""
    variant = VariantNorm(hgvs_c="c.1_2inv", kind="inv", c_start=1, c_end=2, alt=None, meta={})
    with pytest.raises(ValueError, match="Unsupported variant kind"):
        apply_single_variant(ref_seq, variant)