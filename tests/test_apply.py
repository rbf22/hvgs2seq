import pytest
from hgvs2seq.models import VariantNorm, VariantType, GenomicPosition
from hgvs2seq.apply.single import apply_single_variant
from hgvs2seq.apply.batch import apply_variants_in_batch

# A simple reference sequence for testing
# 0-based: 0A 1T 2G 3C 4A 5A 6C 7G 8T 9G 10C 11A 12T
# c. pos:  1A 2T 3G 4C 5A 6A 7C 8G 9T 10G 11C 12A 13T
REF_SEQ = "ATGCAACGTGCAT"

def test_apply_substitution():
    variant = VariantNorm(
        hgvs="NM_123456.7:c.4C>G",
        variant_type=VariantType.SUBSTITUTION,
        transcript_id="NM_123456.7",
        start=4,
        end=4,
        ref="C",
        alt="G",
        hgvs_c="NM_123456.7:c.4C>G"
    )
    expected_seq = "ATGGAACGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_deletion():
    variant = VariantNorm(
        hgvs="NM_123456.7:c.5_6del",
        variant_type=VariantType.DELETION,
        transcript_id="NM_123456.7",
        start=5,
        end=6,
        ref="AA",
        alt="",
        hgvs_c="NM_123456.7:c.5_6del"
    )
    expected_seq = "ATGCCGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_insertion():
    variant = VariantNorm(
        hgvs="NM_123456.7:c.5_6insTT",
        variant_type=VariantType.INSERTION,
        transcript_id="NM_123456.7",
        start=5,
        end=6,
        ref="",
        alt="TT",
        hgvs_c="NM_123456.7:c.5_6insTT"
    )
    expected_seq = "ATGCATTACGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_duplication():
    variant = VariantNorm(
        hgvs="NM_123456.7:c.7_9dup",
        variant_type=VariantType.DUPLICATION,
        transcript_id="NM_123456.7",
        start=7,
        end=9,
        ref="",
        alt="CGT",
        hgvs_c="NM_123456.7:c.7_9dup"
    )
    # Region c.7_9 is CGT. Result should be ATGCAA + CGT + CGT + GCAT
    expected_seq = "ATGCAACGTCGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_inversion():
    variant = VariantNorm(
        hgvs="NM_123456.7:c.2_4inv",
        variant_type=VariantType.INVERSION,
        transcript_id="NM_123456.7",
        start=2,
        end=4,
        ref="TGC",
        alt="GCA",
        hgvs_c="NM_123456.7:c.2_4inv"
    )
    # Region c.2_4 is TGC. rev-comp is GCA. Result should be A + GCA + AACGTGCAT
    expected_seq = "AGCAAACGTGCAT"
    result_seq = apply_single_variant(REF_SEQ, variant)
    assert result_seq == expected_seq

def test_apply_batch_order():
    variants = [
        VariantNorm(
            hgvs="NM_123456.7:c.9T>A",
            variant_type=VariantType.SUBSTITUTION,
            transcript_id="NM_123456.7",
            start=9,
            end=9,
            ref="T",
            alt="A",
            hgvs_c="NM_123456.7:c.9T>A"
        ),
    ]
    # The test is applying c.9T>A (T->A at position 9, 1-based)
    # REF_SEQ:  ATGCAACGTGCAT
    # Position: 123456789 10 11 12 13
    # Change:   T at position 9 to A
    # Expected: ATGCAACGAGCAT
    expected_seq = "ATGCAACGAGCAT"
    result = apply_variants_in_batch(REF_SEQ, variants)
    assert result[0] == expected_seq

def test_apply_batch_phasing():
    variants = [
        VariantNorm(
            hgvs="NM_123456.7:c.2T>A",
            variant_type=VariantType.SUBSTITUTION,
            transcript_id="NM_123456.7",
            start=2,
            end=2,
            ref="T",
            alt="A",
            hgvs_c="NM_123456.7:c.2T>A",
            phase_group=1
        ),
        VariantNorm(
            hgvs="NM_123456.7:c.10G>T",
            variant_type=VariantType.SUBSTITUTION,
            transcript_id="NM_123456.7",
            start=10,
            end=10,
            ref="G",
            alt="T",
            hgvs_c="NM_123456.7:c.10G>T",
            phase_group=1
        ),
        VariantNorm(
            hgvs="NM_123456.7:c.4C>G",
            variant_type=VariantType.SUBSTITUTION,
            transcript_id="NM_123456.7",
            start=4,
            end=4,
            ref="C",
            alt="G",
            hgvs_c="NM_123456.7:c.4C>G",
            phase_group=2
        ),
        VariantNorm(
            hgvs="NM_123456.7:c.9T>A",
            variant_type=VariantType.SUBSTITUTION,
            transcript_id="NM_123456.7",
            start=9,
            end=9,
            ref="T",
            alt="A",
            hgvs_c="NM_123456.7:c.9T>A",
            phase_group=2
        ),
        VariantNorm(
            hgvs="NM_123456.7:c.12A>C",
            variant_type=VariantType.SUBSTITUTION,
            transcript_id="NM_123456.7",
            start=12,
            end=12,
            ref="A",
            alt="C",
            hgvs_c="NM_123456.7:c.12A>C"
        ),
    ]
    expected_hap1 = "AAGCAACGTTCAT"
    expected_hap2 = "ATGGAACGAGCAT"
    expected_hap0 = "ATGCAACGTGCCT"
    result = apply_variants_in_batch(REF_SEQ, variants)
    assert result[0] == expected_hap0
    assert result[1] == expected_hap1
    assert result[2] == expected_hap2

def test_invalid_coordinates():
    variant = VariantNorm(
        hgvs="NM_123456.7:c.20G>T",
        variant_type=VariantType.SUBSTITUTION,
        transcript_id="NM_123456.7",
        start=20,
        end=20,
        ref="G",
        alt="T",
        hgvs_c="NM_123456.7:c.20G>T"
    )
    with pytest.raises(ValueError, match="outside the bounds"):
        apply_single_variant(REF_SEQ, variant)