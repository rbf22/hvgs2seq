from typing import TextIO, List
from ..models import SequenceBundle

def write_fasta(bundle: SequenceBundle, file_handle: TextIO, include: set = {"cdna", "protein"}):
    """
    Writes sequences from a SequenceBundle to a file in FASTA format.

    Args:
        bundle: The SequenceBundle containing the sequences.
        file_handle: A writable file handle.
        include: A set specifying which sequences to include ("cdna", "protein").
    """
    provenance = bundle.provenance
    transcript_id = provenance.get("transcript_id", "unknown_transcript")

    def write_record(seq_id: str, seq: str):
        """Helper to write a single FASTA record."""
        if seq:
            file_handle.write(f">{seq_id}\n")
            # Write sequence, wrapping at 80 characters
            for i in range(0, len(seq), 80):
                file_handle.write(f"{seq[i:i+80]}\n")

    # Write reference sequences
    if "cdna" in include:
        write_record(f"{transcript_id}|ref|cdna", bundle.cdna_ref)
    if "protein" in include and bundle.protein_ref:
        write_record(f"{transcript_id}|ref|protein", bundle.protein_ref)

    # Write edited sequences for each haplotype
    for i, (cdna_edited, protein_edited) in enumerate(zip(bundle.cdna_edited, bundle.protein_edited)):
        haplotype_id = i + 1
        if "cdna" in include:
            write_record(f"{transcript_id}|haplotype={haplotype_id}|cdna", cdna_edited)
        if "protein" in include and protein_edited:
            write_record(f"{transcript_id}|haplotype={haplotype_id}|protein", protein_edited)