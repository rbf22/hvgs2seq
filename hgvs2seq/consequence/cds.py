"""
Handles CDS extraction, translation, and consequence analysis.
"""
import logging
from ..models import TranscriptConfig, ProteinOutcome

_logger = logging.getLogger(__name__)

# Standard DNA codon table
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_sequence(seq: str, stop_symbol='*') -> str:
    """Translates a DNA sequence into a protein sequence, including stop codons."""
    protein = []
    seq_len = len(seq) - (len(seq) % 3)
    for i in range(0, seq_len, 3):
        codon = seq[i:i+3].upper()
        amino_acid = GENETIC_CODE.get(codon, 'X').replace('_', stop_symbol)
        protein.append(amino_acid)
    return "".join(protein)

def analyze_consequences(
    ref_cdna: str,
    edited_cdna: str,
    config: TranscriptConfig
) -> ProteinOutcome:
    """
    Analyzes the consequences of variants on the protein sequence with a clear priority.
    """
    if config.cds_start_c is None or config.cds_end_c is None:
        return ProteinOutcome(protein_sequence=None, hgvs_p=None, consequence="non_coding_transcript")

    cds_start_idx = config.cds_start_c - 1
    cds_end_idx = config.cds_end_c

    ref_cds = ref_cdna[cds_start_idx:cds_end_idx]
    edited_cds = edited_cdna[cds_start_idx:cds_end_idx]

    # Priority 0: Start Loss
    if ref_cds[:3] == 'ATG' and edited_cds[:3] != 'ATG':
        return ProteinOutcome(protein_sequence=translate_sequence(edited_cds), hgvs_p="p.Met1?", consequence="start_loss")

    ref_protein = translate_sequence(ref_cds)
    edited_protein = translate_sequence(edited_cds)

    ref_prot_before_stop = ref_protein.split('*', 1)[0]
    edited_prot_before_stop = edited_protein.split('*', 1)[0]

    # Default values
    consequence = "unknown"
    hgvs_p = "p.?"

    # Priority 1: Stop-gain (most severe outcome)
    if '*' in edited_protein and edited_protein.find('*') < len(ref_prot_before_stop):
        consequence = "stop_gain"
        stop_pos = edited_protein.find('*')
        hgvs_p = f"p.{ref_protein[stop_pos]}{stop_pos+1}*"

    # Priority 2: Frameshift (if not already a stop-gain)
    elif len(ref_cds) % 3 != len(edited_cds) % 3:
        consequence = "frameshift_variant"
        hgvs_p = "p.?" # Simplified for now

    # Priority 3: Stop-loss
    elif ref_protein.endswith('*') and not edited_protein.endswith('*'):
        consequence = "stop_loss"
        hgvs_p = f"p.*{len(ref_prot_before_stop)+1}ext"

    # Priority 4: In-frame indel
    elif len(ref_prot_before_stop) != len(edited_prot_before_stop):
        consequence = "in_frame_indel"
        hgvs_p = "p.?"

    # Priority 5: Missense
    elif ref_prot_before_stop != edited_prot_before_stop:
        consequence = "missense_variant"
        for i, (ref_aa, edit_aa) in enumerate(zip(ref_prot_before_stop, edited_prot_before_stop)):
            if ref_aa != edit_aa:
                hgvs_p = f"p.{ref_aa}{i+1}{edit_aa}"
                break

    # Priority 6: Synonymous or no change
    elif ref_prot_before_stop == edited_prot_before_stop:
        if ref_cds != edited_cds:
            consequence = "synonymous_variant"
        else:
            consequence = "no_change"
        hgvs_p = "p.="

    final_protein_seq = edited_prot_before_stop
    if '*' in edited_protein:
        final_protein_seq += '*'

    return ProteinOutcome(
        protein_sequence=final_protein_seq,
        hgvs_p=hgvs_p,
        consequence=consequence
    )