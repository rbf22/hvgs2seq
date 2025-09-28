import pytest
from hgvs2seq.models import VariantNorm
from hgvs2seq.apply.single import apply_single_variant
from hgvs2seq.apply.batch import apply_variants_in_batch

# A simple reference sequence for testing
# 0-based: 0A 1T 2G 3C 4A 5A 6C 7G 8T 9G 10C 11A 12T
# c. pos:  1A 2T 3G 4C 5A 6A 7C 8G 9T 10G 11C 12A 13T
REF_SEQ = "ATGCAACGTGCAT"

def test_apply_substitution():
    variant = VariantNorm(hgvs_c="c.4C>G", kind="sub", c_start=4, c_end=4, alt="G")
    expected_seq = "ATGGAACGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_deletion():
    variant = VariantNorm(hgvs_c="c.5_6del", kind="del", c_start=5, c_end=6, alt="")
    expected_seq = "ATGCCGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_insertion():
    variant = VariantNorm(hgvs_c="c.5_6insTT", kind="ins", c_start=5, c_end=6, alt="TT")
    expected_seq = "ATGCATTACGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_duplication():
    variant = VariantNorm(hgvs_c="c.7_9dup", kind="dup", c_start=7, c_end=9, alt="")
    # Region c.7_9 is CGT. Result should be ATGCAA + CGT + CGT + GCAT
    expected_seq = "ATGCAACGTCGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_inversion():
    variant = VariantNorm(hgvs_c="c.2_4inv", kind="inv", c_start=2, c_end=4, alt="")
    # Region c.2_4 is TGC. rev-comp is GCA. Result should be A + GCA + ACGTGCAT
    expected_seq = "AGCAACGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_batch_order():
    variants = [
        VariantNorm(hgvs_c="c.4C>G", kind="sub", c_start=4, c_end=4, alt="G"),
        VariantNorm(hgvs_c="c.9T>A", kind="sub", c_start=9, c_end=9, alt="A"),
    ]
    expected_seq = "ATGGAACGAGCAT"
    result = apply_variants_in_batch(REF_SEQ, variants)
    assert result[0] == expected_seq

def test_apply_batch_phasing():
    variants = [
        VariantNorm(hgvs_c="c.2T>A", kind="sub", c_start=2, c_end=2, alt="A", phase_group=1),
        VariantNorm(hgvs_c="c.10G>T", kind="sub", c_start=10, c_end=10, alt="T", phase_group=1),
        VariantNorm(hgvs_c="c.4C>G", kind="sub", c_start=4, c_end=4, alt="G", phase_group=2),
        VariantNorm(hgvs_c="c.9T>A", kind="sub", c_start=9, c_end=9, alt="A", phase_group=2),
        VariantNorm(hgvs_c="c.12A>C", kind="sub", c_start=12, c_end=12, alt="C"),
    ]
    expected_hap1 = "AAGCAACGTTCAT"
    expected_hap2 = "ATGGAACGAGCAT"
    expected_hap0 = "ATGCAACGTGCCT"
    result = apply_variants_in_batch(REF_SEQ, variants)
    assert result[0] == expected_hap0
    assert result[1] == expected_hap1
    assert result[2] == expected_hap2

def test_invalid_coordinates():
    variant = VariantNorm(hgvs_c="c.20G>T", kind="sub", c_start=20, c_end=20, alt="T")
    with pytest.raises(ValueError, match="outside the bounds"):
        apply_single_variant(REF_SEQ, variant)