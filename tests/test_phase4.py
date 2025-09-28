import pytest
from unittest.mock import patch, MagicMock

import hgvs.parser

from hgvs2seq.config import TranscriptConfig
from hgvs2seq.consequence.nmd import check_nmd
from hgvs2seq import apply_variants
from hgvs2seq.parse import VariantIn
from hgvs2seq.refseq import get_reference_sequence


# --- NMD Tests ---

@pytest.fixture
def nmd_test_config() -> TranscriptConfig:
    """A fixture for a transcript with three exons for testing NMD."""
    return TranscriptConfig(
        transcript_id="NM_NMD_TEST",
        gene_symbol="NMDGENE",
        assembly="GRCh38",
        strand=1,
        exons=[(1, 100), (201, 300), (401, 500)],  # 3 exons of 100bp each
        cds_start_c=21,
        cds_end_c=260,  # CDS length is 240 bp (80 codons)
    )


def test_nmd_likely_ptc_far_from_junction(nmd_test_config: TranscriptConfig):
    """Tests a classic NMD case where the PTC is far upstream of the last junction."""
    config = nmd_test_config
    # Last exon-exon junction is at cDNA position 200.
    protein_ref = "A" * 79 + "*"

    # Introduce a PTC at codon 10 (AA position 9).
    # This corresponds to cDNA position 21 + 9*3 = 48.
    # Distance to junction = 200 - 48 = 152. This is >= 50, so NMD is likely.
    edited_cdna = "A" * 20 + ("AAA" * 9) + "TAA" + ("AAA" * 70) + "C" * (300 - 20 - 240)

    result = check_nmd(edited_cdna, protein_ref, config)

    assert result["result"] == "NMD_likely"
    assert result["distance"] == 152


def test_nmd_escaped_ptc_in_last_exon(nmd_test_config: TranscriptConfig):
    """Tests NMD escape when the PTC is in the last exon."""
    config = nmd_test_config
    protein_ref = "A" * 79 + "*"

    # Introduce a PTC in the last exon (i.e., after cDNA position 200).
    # Put it at codon 70 (AA position 69). cDNA pos = 21 + 69*3 = 228.
    edited_cdna = "A" * 20 + ("AAA" * 69) + "TAA" + ("AAA" * 10) + "C" * 40

    result = check_nmd(edited_cdna, protein_ref, config)

    assert result["result"] == "NMD_escaped"
    assert "last exon" in result["reason"]


def test_nmd_escaped_ptc_close_to_junction(nmd_test_config: TranscriptConfig):
    """Tests NMD escape when the PTC is just within the 50-nt threshold."""
    config = nmd_test_config
    protein_ref = "A" * 79 + "*"

    # Last junction at 200. Threshold is 50.
    # Put a PTC at codon 60 (AA pos 59). cDNA pos = 21 + 59*3 = 198.
    # Distance = 200 - 198 = 2. This is < 50, so NMD should be escaped.
    edited_cdna = "A" * 20 + ("AAA" * 59) + "TAA" + ("AAA" * 20) + "C" * 40

    result = check_nmd(edited_cdna, protein_ref, config)

    assert result["result"] == "NMD_escaped"
    assert result["distance"] == 2


# --- SpliceAI Integration Tests ---

@pytest.fixture
def spliceai_test_config() -> TranscriptConfig:
    """A simple config for testing SpliceAI integration."""
    return TranscriptConfig(
        transcript_id="NM_SPLICE_TEST",
        gene_symbol="SPLICEGENE",
        assembly="GRCh38",
        strand=1,
        exons=[(1000, 1100)],
        cds_start_c=1,
        cds_end_c=101,
    )


@patch("hgvs2seq.get_reference_sequence")
@patch("hgvs2seq.splicing.spliceai.project_c_to_g")
@patch("hgvs2seq.splicing.spliceai.subprocess.run")
def test_spliceai_annotation_integration(
    mock_subprocess_run: MagicMock,
    mock_project_c_to_g: MagicMock,
    mock_get_ref_seq: MagicMock,
    spliceai_test_config: TranscriptConfig,
):
    """
    Tests that `apply_variants` correctly calls the SpliceAI annotator
    and integrates the results, mocking the external process and all I/O.
    """
    # 1. Setup Mocks
    hp = hgvs.parser.Parser()
    var_g = hp.parse_hgvs_variant("NC_000001.11:g.1004G>A")
    mock_project_c_to_g.return_value = var_g
    mock_get_ref_seq.return_value = "A" * 4 + "G" + "A" * 96

    def side_effect_run(command: list, *args: list, **kwargs: dict):
        output_vcf_path = command[command.index("-O") + 1]
        with open(output_vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write(
                "NC_000001.11\t1004\tNM_SPLICE_TEST:c.5G>A\tG\tA\t.\t.\t"
                "SpliceAI=A|SPLICEGENE|0.81|0.12|0.03|0.94|2|-10|30|-40\n"
            )
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stderr = "SpliceAI finished successfully."
        return mock_process

    mock_subprocess_run.side_effect = side_effect_run

    # 2. Prepare Inputs
    config = spliceai_test_config
    variants = [VariantIn(hgvs="NM_SPLICE_TEST:c.5G>A")]

    # 3. Call the main function
    bundle = apply_variants(
        cfg=config,
        variants=variants,
        annotate_splicing=True,
        spliceai_params={"reference_path": "/fake/path/to/hg38.fa"},
    )

    # 4. Assertions
    assert not bundle.provenance["edit_plan_warnings"]
    haplotype_annotations = bundle.annotations["haplotypes"][0]
    splicing_results = haplotype_annotations["splicing"]
    splice_scores = splicing_results["NM_SPLICE_TEST:c.5G>A"]

    assert splice_scores["symbol"] == "SPLICEGENE"
    assert splice_scores["ds_ag"] == 0.81
    assert splice_scores["ds_dl"] == 0.94
    assert splice_scores["dp_ag"] == 2
    assert splice_scores["dp_dl"] == -40

    mock_project_c_to_g.assert_called_once()
    mock_subprocess_run.assert_called_once()
    mock_get_ref_seq.assert_called_once_with(config.transcript_id)