import pytest
from hgvs2seq.config import TranscriptConfig
from hgvs2seq.consequence.nmd import check_nmd

# A consistent, well-defined reference sequence for testing.
# Exons (conceptual cDNA locations): 1-100, 101-200, 201-300.
# CDS: c.25-c.279. Length is 255 bp (85 codons).
# The reference CDS has its stop codon at the end.
REF_CDS = ("ATG" * 84) + "TGA"
REF_CDNA = "T" * 24 + REF_CDS + "T" * (300 - 24 - len(REF_CDS))
REF_PROTEIN = "M" * 84 + "*"

@pytest.fixture
def nmd_transcript_config() -> TranscriptConfig:
    """
    Provides a transcript configuration that is consistent with the
    reference sequences defined in this test module.
    """
    return TranscriptConfig(
        transcript_id="NM_NMD",
        gene_symbol="NMD_GENE",
        assembly="GRCh38",
        strand=1,
        # 3 exons, each 100bp long for simplicity in calculating junction positions
        exons=[(1, 100), (201, 300), (401, 500)],
        cds_start_c=25,
        cds_end_c=279,  # 25 + 255 - 1 = 279
    )

def test_nmd_likely_with_default_threshold(nmd_transcript_config: TranscriptConfig):
    """
    Tests if a PTC >50bp from the last junction is correctly flagged as 'NMD_likely'.
    """
    config = nmd_transcript_config
    # Last junction is at cDNA position 200.
    # We will introduce a PTC at c.139 by changing the codon from ATG -> TAG.
    # The cDNA position of the PTC is 139.
    # The distance from the last junction is 200 - 139 = 61 nt.
    # Since 61 >= 50 (default threshold), this should be NMD_likely.
    edited_cdna_list = list(REF_CDNA)
    edited_cdna_list[138] = 'T'  # c.139A>T
    edited_cdna_list[139] = 'A'  # c.140T>A, creating a TAG stop codon
    edited_cdna = "".join(edited_cdna_list)

    result = check_nmd(edited_cdna, REF_PROTEIN, config)

    assert result["result"] == "NMD_likely"
    assert result["ptc_position"] == 139
    assert result["distance"] == 61
    assert ">= 50 nt" in result["reason"]

def test_nmd_escaped_with_custom_threshold(nmd_transcript_config: TranscriptConfig):
    """
    Tests if the same variant can be 'NMD_escaped' by providing a higher custom threshold.
    """
    config = nmd_transcript_config
    # Same setup as before, PTC is 61 nt from the last junction.
    edited_cdna_list = list(REF_CDNA)
    edited_cdna_list[138] = 'T'  # c.139A>T
    edited_cdna_list[139] = 'A'  # c.140T>A, creating a TAG stop codon
    edited_cdna = "".join(edited_cdna_list)

    # Now, we set a custom threshold of 70. Since 61 < 70, it should escape.
    result = check_nmd(edited_cdna, REF_PROTEIN, config, nmd_threshold=70)

    assert result["result"] == "NMD_escaped"
    assert result["distance"] == 61
    assert "< 70 nt" in result["reason"]

def test_nmd_escaped_ptc_in_last_exon(nmd_transcript_config: TranscriptConfig):
    """
    Tests if a PTC in the last exon is correctly flagged as 'NMD_escaped',
    regardless of the threshold.
    """
    config = nmd_transcript_config
    # Last junction is at position 200.
    # Let's introduce a PTC at cDNA position 211 (in the last exon) by changing ATG -> TAG.
    edited_cdna_list = list(REF_CDNA)
    edited_cdna_list[210] = 'T'  # c.211A>T
    edited_cdna_list[211] = 'A'  # c.212T>A
    edited_cdna = "".join(edited_cdna_list)

    # This should escape even with a very small threshold.
    result = check_nmd(edited_cdna, REF_PROTEIN, config, nmd_threshold=10)

    assert result["result"] == "NMD_escaped"
    assert result["reason"] == "PTC is located in the last exon."

def test_nmd_escaped_intronless_transcript():
    """
    Tests if an intronless transcript is correctly flagged as 'NMD_escaped'.
    """
    intronless_config = TranscriptConfig(
        transcript_id="NM_INTLRNLESS",
        gene_symbol="INTLRNLESS_GENE",
        assembly="GRCh38",
        strand=1,
        exons=[(1, 300)],  # Only one exon
        cds_start_c=25,
        cds_end_c=279,
    )
    # This PTC would otherwise be NMD-likely if there were junctions.
    edited_cdna_list = list(REF_CDNA)
    edited_cdna_list[138] = 'T'
    edited_cdna_list[139] = 'A'
    edited_cdna = "".join(edited_cdna_list)

    result = check_nmd(edited_cdna, REF_PROTEIN, intronless_config)

    assert result["result"] == "NMD_escaped"
    assert result["reason"] == "Intronless transcript."