from typing import Optional, Tuple, Dict, Any

from ..config import TranscriptConfig

# Standard genetic code for translation
GENETIC_CODE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}

def extract_cds(cdna: str, config: TranscriptConfig) -> Optional[str]:
    """
    Extracts the coding sequence (CDS) from a cDNA sequence based on the
    transcript configuration.
    """
    if config.cds_start_c is None or config.cds_end_c is None:
        return None
    # Convert from 1-based HGVS coordinates to 0-based Python indices
    start_idx = config.cds_start_c - 1
    end_idx = config.cds_end_c
    return cdna[start_idx:end_idx]

def translate(cds: str) -> str:
    """
    Translates a CDS sequence into a protein sequence using the standard genetic code.
    Translation stops after the first stop codon is encountered. Stops are represented by '*'.
    """
    if not cds:
        return ""
    protein = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i:i+3].upper()
        amino_acid = GENETIC_CODE.get(codon, "X")  # 'X' for unknown codons
        protein.append(amino_acid)
        if amino_acid == "*":
            break  # Stop translation at the first stop codon
    return "".join(protein)

def analyze_consequences(
    protein_ref: Optional[str],
    protein_edited: Optional[str],
) -> Dict[str, Any]:
    """
    Compares reference and edited protein sequences to determine the consequence.
    """
    if protein_ref is None or protein_edited is None:
        return {"consequence": "unknown", "details": "CDS not available or not translated."}

    if protein_ref == protein_edited:
        return {"consequence": "synonymous", "details": "No change in protein sequence."}

    # Check for frameshift by comparing lengths and looking for early termination
    ref_len = len(protein_ref)
    edited_len = len(protein_edited)

    # Find first stop codon
    ref_stop_idx = protein_ref.find("*")
    edited_stop_idx = protein_edited.find("*")

    # Handle lengths considering stop codons
    ref_len_no_stop = ref_len if ref_stop_idx == -1 else ref_stop_idx
    edited_len_no_stop = edited_len if edited_stop_idx == -1 else edited_stop_idx

    if ref_len_no_stop != edited_len_no_stop:
        return {"consequence": "frameshift", "details": f"Protein length changed from {ref_len_no_stop} to {edited_len_no_stop} amino acids."}

    # Check for missense or nonsense
    for i in range(ref_len):
        if i >= edited_len or protein_ref[i] != protein_edited[i]:
            if i < edited_len and protein_edited[i] == "*":
                return {
                    "consequence": "nonsense",
                    "details": f"Stop codon introduced at position {i+1} (p.{protein_ref[i]}{i+1}*)."
                }
            elif i >= edited_len:
                 return {
                    "consequence": "frameshift", # Should have been caught by length check, but as a fallback
                    "details": "Protein is truncated."
                }
            else:
                return {
                    "consequence": "missense",
                    "details": f"Amino acid change at position {i+1} from {protein_ref[i]} to {protein_edited[i]} (p.{protein_ref[i]}{i+1}{protein_edited[i]})."
                }

    return {"consequence": "no_protein_change_detected", "details": "Sequences differ but no specific consequence found by this basic analysis."}


def get_protein_sequence_and_consequences(
    cdna_ref: str,
    cdna_edited: str,
    config: TranscriptConfig,
) -> Tuple[Optional[str], Optional[str], Dict[str, Any]]:
    """
    Main function to extract CDS, translate, and analyze consequences.
    """
    cds_ref = extract_cds(cdna_ref, config)
    cds_edited = extract_cds(cdna_edited, config)

    if cds_ref is None or cds_edited is None:
        return None, None, {"consequence": "unknown", "details": "Could not extract CDS."}

    protein_ref = translate(cds_ref)
    protein_edited = translate(cds_edited)

    consequences = analyze_consequences(protein_ref, protein_edited)

    return protein_ref, protein_edited, consequences