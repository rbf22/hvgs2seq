import pytest
from hgvs2seq.config import TranscriptConfig
from hgvs2seq.parse import VariantIn, normalize_variant, VariantNorm

# Using BRCA1 transcript as a stable, well-known example for testing projection.
# The exon/CDS details are simplified for testing but are based on the real transcript.
@pytest.fixture
def brca1_config():
    """Provides a sample TranscriptConfig for BRCA1 (NM_007294.3)."""
    return TranscriptConfig(
        transcript_id="NM_007294.3",
        gene_symbol="BRCA1",
        assembly="GRCh38",
        strand=-1,
        # Exon coordinates are not strictly needed for c. projection but are good practice to have.
        exons=[(43044295, 43045813), (43047631, 43047723)], # Simplified for this example
        cds_start_c=225, # Based on real data
        cds_end_c=5792,  # Based on real data
    )

def test_normalize_genomic_snv(brca1_config):
    """Tests projecting a simple genomic SNV to cDNA coordinates."""
    # This variant corresponds to a known pathogenic BRCA1 variant.
    # The reference base at this position is A, so the HGVS string must reflect that.
    variant_in = VariantIn(hgvs="NC_000017.11:g.43091523A>T")

    norm_variant = normalize_variant(variant_in, brca1_config)

    assert isinstance(norm_variant, VariantNorm)
    # The projection to the reverse strand transcript correctly inverts the alleles.
    # The expected position is c.4008, as determined by the UTA database.
    assert norm_variant.hgvs_c == "NM_007294.3:c.4008T>A"
    assert norm_variant.kind == "sub"
    assert norm_variant.c_start == 4008
    assert norm_variant.c_end == 4008
    assert norm_variant.alt == "A"
    assert norm_variant.meta["original_hgvs"] == variant_in.hgvs

def test_normalize_cDNA_deletion(brca1_config):
    """Tests a variant already in cDNA coordinates (a small deletion)."""
    # A well-known 2-base deletion in BRCA1
    hgvs_string = "NM_007294.3:c.68_69del"
    variant_in = VariantIn(hgvs=hgvs_string)

    norm_variant = normalize_variant(variant_in, brca1_config)

    assert norm_variant.hgvs_c == hgvs_string
    assert norm_variant.kind == "del"
    assert norm_variant.c_start == 68
    assert norm_variant.c_end == 69
    assert norm_variant.alt is None

def test_normalize_invalid_hgvs_string(brca1_config):
    """Tests that a malformed HGVS string raises a ValueError."""
    variant_in = VariantIn(hgvs="this is not a valid hgvs string")

    with pytest.raises(ValueError, match="Failed to parse HGVS string"):
        normalize_variant(variant_in, brca1_config)

def test_normalize_unprojectable_variant(brca1_config):
    """Tests that a variant on a different chromosome fails projection."""
    # A variant on chromosome 1 cannot be projected to BRCA1 on chr 17.
    variant_in = VariantIn(hgvs="NC_000001.11:g.12345A>T")

    with pytest.raises(ValueError, match="Failed to project variant"):
        normalize_variant(variant_in, brca1_config)

def test_normalize_cDNA_insertion(brca1_config):
    """Tests a cDNA insertion."""
    hgvs_string = "NM_007294.3:c.5258_5259insA"
    variant_in = VariantIn(hgvs=hgvs_string)

    norm_variant = normalize_variant(variant_in, brca1_config)

    assert norm_variant.hgvs_c == hgvs_string
    assert norm_variant.kind == "ins"
    assert norm_variant.c_start == 5258
    assert norm_variant.c_end == 5259
    assert norm_variant.alt == "A"