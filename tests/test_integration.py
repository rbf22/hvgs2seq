import pytest
import json
from hgvs2seq.config import TranscriptConfig
from hgvs2seq.parse import VariantNorm
from hgvs2seq.models import EditPlan, SequenceBundle
from hgvs2seq.apply.batch import apply_edit_plan
from hgvs2seq.consequence.cds import get_protein_sequence_and_consequences
from hgvs2seq.io.fasta import write_fasta
from hgvs2seq.io.jsonio import write_json
import io

# A sample transcript configuration for a toy gene
@pytest.fixture
def sample_transcript_config() -> TranscriptConfig:
    return TranscriptConfig(
        transcript_id="NM_001",
        gene_symbol="TOYGENE",
        assembly="GRCh38",
        strand=1,
        exons=[(101, 200), (301, 400)],
        cds_start_c=10,  # CDS starts at the 10th base of the cDNA
        cds_end_c=109,   # CDS ends at the 109th base
    )

# Sample reference cDNA sequence corresponding to the config
# Let's make the CDS region easy to track:
# cDNA:   AAAAAAAAG|ATG GGC AAT GTC GAG TGA|AAAAAAAAAA...
# 1-based: 1       9|10                 109|110
# Protein:         |M  G  N  V  E  *  |
REF_CDNA = "AAAAAAAAG" + "ATGGGCAATGTCGAGTGA" + "AAAAAAAAA" * 10

def test_end_to_end_missense_variant(sample_transcript_config: TranscriptConfig):
    """
    Tests the full pipeline with a single missense variant.
    It verifies sequence editing, CDS extraction, translation, consequence analysis, and output generation.
    """
    config = sample_transcript_config

    # This variant will change c.13G>C (in the CDS), which is GGC -> CGC, so Gly -> Arg (G -> R)
    missense_variant = VariantNorm(
        hgvs_c="NM_001:c.13G>C",
        kind="sub",
        c_start=13,
        c_end=13,
        alt="C",
        meta={"original_hgvs": "NM_001:c.13G>C"}
    )

    # 1. Create an EditPlan
    edit_plan = EditPlan(
        haplotypes=[[missense_variant]],
        policy="order_by_pos",
        warnings=[]
    )

    # 2. Apply the edit plan
    edited_cdn_sequences = apply_edit_plan(REF_CDNA, edit_plan)
    assert len(edited_cdn_sequences) == 1
    edited_cdna = edited_cdn_sequences[0]

    # Expected cDNA: The 'G' at index 12 (1-based 13) should become 'C'
    # Original CDS part: ATGGGCAAT...
    # Edited CDS part:   ATGCGCAAT...
    expected_cdna = "AAAAAAAAG" + "ATGCGCAATGTCGAGTGA" + "AAAAAAAAA" * 10
    assert edited_cdna == expected_cdna

    # 3. Get protein sequence and consequences
    protein_ref, protein_edited, consequences = get_protein_sequence_and_consequences(
        REF_CDNA, edited_cdna, config
    )

    assert protein_ref == "MGNVE*"
    assert protein_edited == "MRNVE*"
    assert consequences["consequence"] == "missense"
    assert "p.G2R" in consequences["details"]

    # 4. Create a SequenceBundle
    bundle = SequenceBundle(
        cdna_ref=REF_CDNA,
        cdna_edited=[edited_cdna],
        mrna_ref=REF_CDNA, # Assuming same for this test
        mrna_edited=[edited_cdna],
        protein_ref=protein_ref,
        protein_edited=[protein_edited],
        annotations={"haplotype_1": consequences},
        provenance={"transcript_id": config.transcript_id, "variants_applied": [missense_variant.hgvs_c]}
    )

    # 5. Test output generation
    # Test FASTA output
    fasta_io = io.StringIO()
    write_fasta(bundle, fasta_io)
    fasta_output = fasta_io.getvalue()
    assert ">NM_001|ref|protein" in fasta_output
    assert "MGNVE*" in fasta_output
    assert ">NM_001|haplotype=1|protein" in fasta_output
    assert "MRNVE*" in fasta_output

    # Test JSON output
    json_io = io.StringIO()
    write_json(bundle, json_io)
    json_output = json_io.getvalue()
    bundle_dict = json.loads(json_output)
    assert bundle_dict["protein_edited"][0] == "MRNVE*"
    assert bundle_dict["annotations"]["haplotype_1"]["consequence"] == "missense"