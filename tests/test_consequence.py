import pytest
from hgvs2seq.models import TranscriptConfig, ProteinOutcome
from hgvs2seq.consequence.cds import analyze_consequences
from hgvs2seq.consequence.nmd import check_nmd

@pytest.fixture
def transcript_config():
    from hgvs2seq.config import load_config
    return load_config("tests/fixtures/test_transcript_config.json")

# A reference cDNA that is consistent with the transcript config fixture.
# 5' UTR (10bp) + CDS (186bp) + 3' UTR (19bp) = 215bp total
REF_CDNA = "GATTACAGAT" + "ATG" + "GGC" * 30 + "TCC" * 30 + "TAA" + "GATTACAGATTACAGATTA"

def get_edited_cdna(cds_change_func):
    """Helper to create an edited cDNA from a function that modifies the CDS."""
    ref_cds = REF_CDNA[10:196]
    edited_cds = cds_change_func(ref_cds)
    return REF_CDNA[:10] + edited_cds + REF_CDNA[196:]

def test_missense(transcript_config):
    # GGC (Gly) -> GAC (Asp)
    edited_cdna = get_edited_cdna(lambda cds: cds[:3] + "GAC" + cds[6:])
    outcome = analyze_consequences(REF_CDNA, edited_cdna, transcript_config)
    assert outcome.consequence == "missense_variant"

def test_nonsense(transcript_config):
    # GGC (Gly) -> TGA (Stop)
    edited_cdna = get_edited_cdna(lambda cds: cds[:3] + "TGA" + cds[6:])
    outcome = analyze_consequences(REF_CDNA, edited_cdna, transcript_config)
    assert outcome.consequence == "stop_gain"

def test_frameshift_causes_stop_gain(transcript_config):
    # Delete one base, causing a frameshift that results in an early stop.
    # The most severe consequence is stop_gain.
    edited_cdna = get_edited_cdna(lambda cds: cds[:4] + cds[5:])
    outcome = analyze_consequences(REF_CDNA, edited_cdna, transcript_config)
    assert outcome.consequence == "stop_gain"

def test_in_frame_deletion(transcript_config):
    # Delete a whole codon.
    edited_cdna = get_edited_cdna(lambda cds: cds[:3] + cds[6:])
    outcome = analyze_consequences(REF_CDNA, edited_cdna, transcript_config)
    assert outcome.consequence == "in_frame_indel"

def test_nmd_likely(transcript_config):
    # Introduce a stop codon early in the sequence.
    edited_cdna = get_edited_cdna(lambda cds: cds[:3] + "TAA" + cds[6:])
    protein_outcome = analyze_consequences(REF_CDNA, edited_cdna, transcript_config)
    assert protein_outcome.consequence == "stop_gain"
    nmd_outcome = check_nmd(protein_outcome, transcript_config)
    assert nmd_outcome.status == "likely"

def test_nmd_escape(transcript_config):
    # Introduce a stop codon in the last exon.
    edited_cdna = get_edited_cdna(lambda cds: cds[:120] + "TGA" + cds[123:])
    protein_outcome = analyze_consequences(REF_CDNA, edited_cdna, transcript_config)
    assert protein_outcome.consequence == "stop_gain"
    nmd_outcome = check_nmd(protein_outcome, transcript_config)
    assert nmd_outcome.status == "escape"