"""
The command-line interface for the hgvs2seq tool.
"""
import click
import logging
import sys
from typing import List, Optional
import re

from . import __version__
from .models import VariantIn, SequenceBundle, TranscriptOutcome
from .config import load_config
from .refseq import get_reference_cDNA
from .parse import parse_and_normalize_variants
from .apply.batch import apply_variants_in_batch
from .consequence.cds import analyze_consequences
from .consequence.nmd import check_nmd
from .splicing.spliceai import analyze_all_variants_for_splicing
from .io.jsonio import generate_json_output
from .io.fasta import generate_fasta_output

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
_logger = logging.getLogger(__name__)

def parse_variants_file(file_path: str) -> List[VariantIn]:
    """Parses a file of HGVS strings, one per line, with optional phase comments."""
    variants = []
    phase_regex = re.compile(r'#\s*phase\s*=\s*(\d+)')
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            phase_group: Optional[int] = None
            match = phase_regex.search(line)
            if match:
                phase_group = int(match.group(1))
                # Remove comment from hgvs string
                line = phase_regex.sub('', line).strip()

            variants.append(VariantIn(hgvs=line, phase_group=phase_group))
    return variants

@click.command()
@click.version_option(version=__version__)
@click.option('--config-path', '-c', required=True, type=click.Path(exists=True, dir_okay=False),
              help='Path to the transcript configuration JSON file.')
@click.option('--variants-path', '-v', required=True, type=click.Path(exists=True, dir_okay=False),
              help='Path to a text file with one HGVS variant per line.')
@click.option('--json-out', '-jo', type=click.Path(dir_okay=False),
              help='Path to write the output JSON file.')
@click.option('--fasta-out', '-fo', type=click.Path(dir_okay=False),
              help='Path to write the output FASTA file.')
@click.option('--emit-protein/--no-emit-protein', default=True, help='Include protein sequences in FASTA output.')
@click.option('--emit-cdna/--no-emit-cdna', default=True, help='Include cDNA sequences in FASTA output.')
def main(config_path, variants_path, json_out, fasta_out, emit_protein, emit_cdna):
    """
    Takes a list of HGVS variants and a transcript, and outputs the resulting sequences.
    """
    _logger.info(f"hgvs2seq v{__version__} starting...")

    try:
        # 1. Load configuration and reference data
        config = load_config(config_path)
        variants_in = parse_variants_file(variants_path)
        ref_cdna = get_reference_cDNA(config.transcript_id)

        # 2. Parse and normalize variants
        norm_variants = parse_and_normalize_variants(variants_in, config)

        # 3. Apply variants in batch
        edited_sequences_by_haplotype = apply_variants_in_batch(ref_cdna, norm_variants)

        # 4. Analyze consequences for each haplotype
        primary_outcomes = []
        for hap_id, edited_cdna in edited_sequences_by_haplotype.items():
            _logger.info(f"Analyzing consequences for haplotype {hap_id}...")
            protein_outcome = analyze_consequences(ref_cdna, edited_cdna, config)
            nmd_outcome = check_nmd(protein_outcome, config)

            transcript_outcome = TranscriptOutcome(
                haplotype_id=hap_id,
                mrna_sequence=edited_cdna,
                protein_outcome=protein_outcome,
                nmd=nmd_outcome,
            )
            primary_outcomes.append(transcript_outcome)

        # 5. Perform splicing analysis (annotation)
        splice_annotations = analyze_all_variants_for_splicing(norm_variants, config)

        # 6. Assemble the final SequenceBundle
        bundle = SequenceBundle(
            primary_outcomes=primary_outcomes,
            provenance={
                "tool_version": __version__,
                "transcript_id": config.transcript_id,
                "assembly": config.assembly,
                "input_variants": [v.hgvs for v in variants_in],
            },
            warnings=splice_annotations
        )

        # 7. Generate and write outputs
        if json_out:
            json_output = generate_json_output(bundle)
            with open(json_out, 'w') as f:
                f.write(json_output)
            _logger.info(f"Wrote JSON output to {json_out}")

        if fasta_out:
            fasta_output = generate_fasta_output(bundle, config, emit_protein, emit_cdna)
            with open(fasta_out, 'w') as f:
                f.write(fasta_output)
            _logger.info(f"Wrote FASTA output to {fasta_out}")

        _logger.info("hgvs2seq finished successfully.")

    except Exception as e:
        _logger.error(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == '__main__':
    main()