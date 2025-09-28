"""
Handles generation of FASTA formatted outputs.
"""
import logging
from typing import List, Union
from textwrap import wrap
from ..models import SequenceBundle, TranscriptOutcome, TranscriptConfig

_logger = logging.getLogger(__name__)

def _format_fasta_entry(header: str, sequence: str, line_length: int = 60) -> str:
    """Formats a single sequence into a FASTA entry."""
    wrapped_sequence = "\n".join(wrap(sequence, line_length))
    return f">{header}\n{wrapped_sequence}\n"

def generate_fasta_output(
    bundle: SequenceBundle,
    config: TranscriptConfig,
    emit_protein: bool = True,
    emit_cdna: bool = True
) -> str:
    """
    Generates a FASTA-formatted string containing the edited sequences.

    Args:
        bundle: The SequenceBundle containing the outcomes.
        config: The TranscriptConfig for metadata.
        emit_protein: Whether to include protein sequences in the output.
        emit_cdna: Whether to include cDNA (mRNA) sequences in the output.

    Returns:
        A string in FASTA format.
    """
    fasta_entries = []

    all_outcomes = bundle.primary_outcomes + bundle.alternate_outcomes

    for outcome in all_outcomes:
        # Common part of the header
        base_header = (
            f"{config.transcript_id}|haplotype={outcome.haplotype_id}|"
            f"scenario={outcome.scenario_id}"
        )

        # Add cDNA sequence if requested
        if emit_cdna:
            cdna_header = f"{base_header}|type=cDNA"
            fasta_entries.append(_format_fasta_entry(cdna_header, outcome.mrna_sequence))

        # Add protein sequence if requested and available
        if emit_protein and outcome.protein_outcome and outcome.protein_outcome.protein_sequence:
            protein_header = f"{base_header}|type=protein|consequence={outcome.protein_outcome.consequence}"
            fasta_entries.append(
                _format_fasta_entry(protein_header, outcome.protein_outcome.protein_sequence)
            )

    _logger.info(f"Generated FASTA output with {len(fasta_entries)} entries.")
    return "".join(fasta_entries)