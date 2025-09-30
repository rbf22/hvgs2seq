"""Integration tests for the hgvs2seq package."""
import pytest
from hgvs2seq.models import VariantNorm, TranscriptConfig, VariantType, GenomicPosition
from hgvs2seq.project import project_variant, project_c_to_g

@pytest.fixture
def test_transcript():
    """Return a test transcript configuration."""
    return TranscriptConfig(
        transcript_id="NM_000000.0",
        gene_symbol="TEST",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=1000,
        tx_end=2000,
        exons=[(1000, 1200), (1500, 1700), (1800, 2000)],
        cds_start=1100,
        cds_end=1900
    )

def test_end_to_end_variant_projection(test_transcript):
    """Test the full workflow from genomic to cDNA to protein."""
    # Create a genomic variant
    genomic_variant = VariantNorm(
        hgvs="chr1:g.1050A>G",
        variant_type=VariantType.SUBSTITUTION,
        transcript_id=test_transcript.transcript_id,
        start=1050,
        end=1050,
        ref="A",
        alt="G",
        hgvs_c="",
        hgvs_p=""
    )
    
    # Project to cDNA
    cdna_variant = project_variant(genomic_variant, test_transcript)
    assert cdna_variant.hgvs_c.startswith('c.')
    assert cdna_variant.transcript_id == test_transcript.transcript_id
    
    # Project back to genomic
    genomic_result = project_c_to_g(cdna_variant, test_transcript)
    assert genomic_result.genomic_pos is not None
    assert genomic_result.hgvs.startswith('chr1:g.')
    assert 'A>G' in genomic_result.hgvs

def test_indel_variant(test_transcript):
    """Test handling of insertion/deletion variants."""
    # Test deletion
    del_variant = VariantNorm(
        hgvs="chr1:g.1050_1052del",
        variant_type=VariantType.DELETION,
        transcript_id=test_transcript.transcript_id,
        start=1050,
        end=1052,
        ref="ATG",
        alt="",
        hgvs_c="",
        hgvs_p=""
    )
    
    cdna_del = project_variant(del_variant, test_transcript)
    assert "del" in cdna_del.hgvs_c.lower()
    
    # Test insertion
    ins_variant = VariantNorm(
        hgvs="chr1:g.1050_1051insT",
        variant_type=VariantType.INSERTION,
        transcript_id=test_transcript.transcript_id,
        start=1050,
        end=1051,
        ref="",
        alt="T",
        hgvs_c="",
        hgvs_p=""
    )
    
    cdna_ins = project_variant(ins_variant, test_transcript)
    assert "ins" in cdna_ins.hgvs_c.lower()

def test_exon_boundary_cases(test_transcript):
    """Test variants at exon boundaries."""
    # Variant at the start of an exon
    exon_start_variant = VariantNorm(
        hgvs=f"chr1:g.{test_transcript.exons[0][0]}A>G",
        variant_type=VariantType.SUBSTITUTION,
        transcript_id=test_transcript.transcript_id,
        start=test_transcript.exons[0][0],
        end=test_transcript.exons[0][0],
        ref="A",
        alt="G",
        hgvs_c="",
        hgvs_p=""
    )
    
    cdna_variant = project_variant(exon_start_variant, test_transcript)
    assert cdna_variant.hgvs_c.startswith('c.')
