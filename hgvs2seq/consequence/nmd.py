from typing import Dict, Any, Optional

from ..config import TranscriptConfig
from .cds import extract_cds, translate

# The default distance from a PTC to the final exon-exon junction that typically triggers NMD.
NMD_THRESHOLD_DEFAULT = 50

def check_nmd(
    edited_cdna: str,
    protein_ref: Optional[str],
    config: TranscriptConfig,
    *,
    nmd_threshold: int = NMD_THRESHOLD_DEFAULT,
) -> Dict[str, Any]:
    """
    Applies a standard set of rules to determine if a transcript with a
    premature termination codon (PTC) is a likely candidate for
    Nonsense-Mediated Decay (NMD).

    The primary rule is:
    - If a PTC is located more than 50-55 nucleotides upstream of the
      final exon-exon junction, the transcript is likely targeted for NMD.

    Exceptions:
    - Transcripts with only one exon (intronless).
    - PTCs located in the last exon.

    Args:
        edited_cdna: The full cDNA sequence of the edited transcript.
        protein_ref: The reference protein sequence, including the stop codon (*).
        config: The transcript configuration, containing exon structure and CDS info.

    Returns:
        A dictionary containing the NMD prediction ("NMD_likely", "NMD_escaped", etc.)
        and a rationale for the decision.
    """
    # Exception 1: Intronless transcripts are not subject to NMD by this rule.
    if len(config.exons) <= 1:
        return {"result": "NMD_escaped", "reason": "Intronless transcript."}

    # Determine the position of the last exon-exon junction in cDNA coordinates.
    # This assumes exons in config are ordered by transcription.
    try:
        exon_lengths = [end - start + 1 for start, end in config.exons]
        last_junction_pos_c = sum(exon_lengths[:-1])
    except (TypeError, IndexError):
        return {"result": "unknown", "reason": "Could not determine exon structure from config."}

    # Extract the edited CDS and translate it to find the PTC.
    cds_edited = extract_cds(edited_cdna, config)
    if not cds_edited:
        return {"result": "unknown", "reason": "Could not extract CDS from edited cDNA."}

    protein_edited = translate(cds_edited)
    ptc_aa_pos = protein_edited.find('*')

    # If there's no stop codon at all, NMD is not applicable.
    if ptc_aa_pos == -1:
        return {"result": "NMD_not_applicable", "reason": "No stop codon found in edited sequence."}

    # Check if the stop codon is premature by comparing to the reference.
    if protein_ref is None:
        return {"result": "unknown", "reason": "Reference protein not provided to check for PTC."}

    ref_stop_aa_pos = protein_ref.find('*')
    is_ptc = (ref_stop_aa_pos == -1) or (ptc_aa_pos < ref_stop_aa_pos)

    if not is_ptc:
        return {"result": "NMD_not_applicable", "reason": "Not a premature termination codon (PTC)."}

    # Calculate the 1-based cDNA position of the PTC.
    if config.cds_start_c is None:
        return {"result": "unknown", "reason": "CDS start position not defined in config."}
    ptc_cdna_pos_1based = config.cds_start_c + (ptc_aa_pos * 3)

    # Exception 2: PTC is in the last exon.
    if ptc_cdna_pos_1based > last_junction_pos_c:
        return {
            "result": "NMD_escaped",
            "reason": "PTC is located in the last exon.",
            "ptc_position": ptc_cdna_pos_1based,
            "last_junction_position": last_junction_pos_c,
        }

    # Main NMD rule: Check distance from PTC to the last junction.
    distance = last_junction_pos_c - ptc_cdna_pos_1based

    if distance >= nmd_threshold:
        return {
            "result": "NMD_likely",
            "reason": f"PTC is {distance} nt upstream of the final exon-exon junction (>= {nmd_threshold} nt).",
            "ptc_position": ptc_cdna_pos_1based,
            "last_junction_position": last_junction_pos_c,
            "distance": distance,
        }
    else:
        return {
            "result": "NMD_escaped",
            "reason": f"PTC is {distance} nt upstream of the final exon-exon junction (< {nmd_threshold} nt).",
            "ptc_position": ptc_cdna_pos_1based,
            "last_junction_position": last_junction_pos_c,
            "distance": distance,
        }