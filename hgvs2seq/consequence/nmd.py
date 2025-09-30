"""
Handles prediction of nonsense-mediated decay (NMD).
"""
import logging
from ..models import TranscriptConfig, ProteinOutcome, NMDOutcome

_logger = logging.getLogger(__name__)

# The distance from a PTC to the last exon-exon junction to trigger NMD.
NMD_RULE_DISTANCE_NT = 55

def check_nmd(
    protein_outcome: ProteinOutcome,
    config: TranscriptConfig
) -> NMDOutcome:
    """
    Checks if a transcript with a premature termination codon (PTC) is likely
    to be targeted by nonsense-mediated decay (NMD).

    Args:
        protein_outcome: The outcome of the protein analysis.
        config: The transcript configuration containing exon information.

    Returns:
        An NMDOutcome object with the status and rationale.
    """
    # 1. NMD is only relevant if there's a premature stop codon (stop_gain).
    if protein_outcome.consequence != "stop_gain":
        return NMDOutcome(status="not_applicable", rationale="No premature termination codon found.")

    # 2. Handle intronless transcripts (single exon).
    if len(config.exons) <= 1:
        return NMDOutcome(status="escape", rationale="Transcript is intronless.")

    # 3. Find the position of the PTC in the cDNA.
    # The protein sequence includes the stop codon ('*').
    ptc_protein_pos = protein_outcome.protein_sequence.find('*')
    if ptc_protein_pos == -1:
        # This case should not happen if consequence is stop_gain, but as a safeguard:
        return NMDOutcome(status="not_applicable", rationale="Could not locate PTC in protein sequence.")

    # Convert protein position (0-based) to cDNA position (1-based).
    # The start of the codon is at protein_pos * 3.
    # The position of the stop codon is relative to the CDS start.
    if config.cds_start is None:
        return NMDOutcome(status="not_applicable", rationale="CDS start is not defined.")

    cds_start_idx = config.cds_start - 1
    ptc_cds_pos = ptc_protein_pos * 3
    ptc_cdna_pos = cds_start_idx + ptc_cds_pos + 1 # +1 for 1-based coordinate

    # 4. Calculate the position of the last exon-exon junction in cDNA coordinates.
    # The junction position is the length of all exons except the last one.
    exon_lengths = [(end - start + 1) for start, end in config.exons]
    last_junction_cdna_pos = sum(exon_lengths[:-1])

    # 5. Apply the NMD rule.
    # The rule applies if the PTC is *upstream* of the junction.
    if ptc_cdna_pos > last_junction_cdna_pos:
        # PTC is in the last exon
        status = "escape"
        rationale = (f"PTC at c.{ptc_cdna_pos} is located in the last exon, "
                     f"downstream of the final exon-exon junction at c.{last_junction_cdna_pos}.")
    else:
        distance_to_last_junction = last_junction_cdna_pos - ptc_cdna_pos
        if distance_to_last_junction >= NMD_RULE_DISTANCE_NT:
            status = "likely"
            rationale = (f"PTC at c.{ptc_cdna_pos} is {distance_to_last_junction} nt "
                         f"upstream of the last exon-exon junction (at c.{last_junction_cdna_pos}), "
                         f"which is >= {NMD_RULE_DISTANCE_NT} nt.")
        else:
            status = "escape"
            rationale = (f"PTC at c.{ptc_cdna_pos} is {distance_to_last_junction} nt "
                         f"upstream of the last exon-exon junction (at c.{last_junction_cdna_pos}), "
                         f"which is < {NMD_RULE_DISTANCE_NT} nt.")

    return NMDOutcome(status=status, rationale=rationale)