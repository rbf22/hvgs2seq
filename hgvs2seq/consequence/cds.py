"""
Handles CDS extraction, translation, and consequence analysis.
"""
import logging
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any

@dataclass
class ConsequenceCheckResult:
    """Result of a consequence check."""
    consequence: Optional[str] = None
    hgvs_p: Optional[str] = None
    protein_sequence: Optional[str] = None
    details: Dict[str, Any] = None
    
    def is_positive(self) -> bool:
        """Return True if this check found a consequence."""
        return self.consequence is not None

from ..models import TranscriptConfig, ProteinOutcome

# Standard stop codons
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

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
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_sequence(seq: str, stop_symbol: str = '*') -> str:
    """
    Translate a DNA sequence into a protein sequence.
    
    Args:
        seq: DNA sequence to translate
        stop_symbol: Symbol to use for stop codons (default: '*')
        
    Returns:
        Translated protein sequence with stop symbols
    """
    _logger = logging.getLogger(__name__)
    protein = []
    
    # Log input sequence for debugging
    _logger.debug(f"Translating sequence (first 100bp): {seq[:100] if seq else 'None'}")
    _logger.debug(f"Sequence length: {len(seq) if seq else 0}")
    
    # Check for empty or invalid sequence
    if not seq:
        _logger.error("Empty sequence provided for translation")
        return ""
    
    # Convert to uppercase and validate characters
    seq = seq.upper()
    valid_bases = set('ACGTN')
    invalid_bases = set(seq) - valid_bases
    
    if invalid_bases:
        _logger.warning(f"Invalid bases in sequence: {invalid_bases}")
    
    # Translate each codon
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        
        # Handle incomplete codons at the end of the sequence
        if len(codon) < 3:
            _logger.warning(f"Incomplete codon at position {i}: {codon}")
            break
        
        # Look up the amino acid in the genetic code
        aa = GENETIC_CODE.get(codon, 'X')
        
        # Replace stop symbol if needed
        if aa == '*' and stop_symbol != '*':
            aa = stop_symbol
            
        protein.append(aa)
        
        # Log the first few translations for debugging
        if i < 30:  # Log first 10 codons
            _logger.debug(f"Translated {codon} -> {aa}")
    
    protein_seq = ''.join(protein)
    
    # Log the results
    _logger.debug(f"Translated protein sequence (first 100aa): {protein_seq[:100]}")
    _logger.debug(f"Protein sequence length: {len(protein_seq)}")
    
    return protein_seq

def check_start_loss(
    ref_cds: str,
    edited_cds: str,
    ref_protein: str,
    edited_protein: str
) -> ConsequenceCheckResult:
    """Check for start codon loss."""
    result = ConsequenceCheckResult()
    
    # Check for start codon loss (must be first check)
    if ref_cds.startswith('ATG') and not edited_cds.startswith('ATG'):
        result.consequence = "start_loss"
        result.hgvs_p = "p.Met1?"
        result.protein_sequence = edited_protein
        
    return result

def check_in_frame_indels(
    ref_cds: str,
    edited_cds: str,
    ref_protein: str,
    edited_protein: str
) -> ConsequenceCheckResult:
    """Check for in-frame insertions and deletions."""
    result = ConsequenceCheckResult()
    is_indel = len(ref_cds) != len(edited_cds)
    
    # Only check for in-frame indels if the lengths differ by a multiple of 3
    if not is_indel or (len(ref_cds) - len(edited_cds)) % 3 != 0:
        return result
    
    # Find the position where the sequences first differ
    diff_pos = 0
    for i in range(min(len(ref_cds), len(edited_cds))):
        if i >= len(ref_cds) or i >= len(edited_cds) or ref_cds[i] != edited_cds[i]:
            diff_pos = i
            break
    
    # Get the codon position (1-based)
    aa_pos = diff_pos // 3 + 1
    
    # For in-frame deletions
    if len(edited_cds) < len(ref_cds):
        del_len = len(ref_cds) - len(edited_cds)
        deleted_nt = ref_cds[diff_pos:diff_pos + del_len]
        deleted_aa = translate_sequence(deleted_nt)
        
        result.consequence = "in_frame_indel"
        result.hgvs_p = f"p.{deleted_aa[0] if deleted_aa else '?'}{aa_pos}del"
        result.protein_sequence = edited_protein
    
    # For in-frame insertions
    elif len(edited_cds) > len(ref_cds):
        ins_pos = diff_pos
        inserted_nt = edited_cds[ins_pos:ins_pos + (len(edited_cds) - len(ref_cds))]
        inserted_aa = translate_sequence(inserted_nt)
        
        result.consequence = "in_frame_indel"
        result.hgvs_p = f"p.{ref_protein[aa_pos-1] if aa_pos <= len(ref_protein) else '?'}{aa_pos}ins{inserted_aa}"
        result.protein_sequence = edited_protein
    
    return result

def check_frameshifts(
    ref_cds: str,
    edited_cds: str,
    ref_protein: str,
    edited_protein: str
) -> ConsequenceCheckResult:
    """Check for frameshift variants."""
    _logger = logging.getLogger(__name__)
    result = ConsequenceCheckResult()
    
    # Check if the CDS sequences are different lengths (indel)
    is_indel = len(ref_cds) != len(edited_cds)
    
    # Calculate the length difference (positive for deletion, negative for insertion)
    length_diff = len(ref_cds) - len(edited_cds)
    
    # This is a frameshift if:
    # 1. It's an indel and the length difference is not a multiple of 3, or
    # 2. The protein sequences have different lengths (excluding terminal stop codons)
    ref_protein_trimmed = ref_protein.rstrip('*')
    edited_protein_trimmed = edited_protein.rstrip('*')
    
    # Check if this might be a stop_lost variant (handled by check_stop_gain)
    ref_has_stop = '*' in ref_protein or 'X' in ref_protein
    edited_has_stop = '*' in edited_protein or 'X' in edited_protein
    
    # If it's an indel and not a multiple of 3, it's definitely a frameshift
    if is_indel and abs(length_diff) % 3 != 0:
        is_frameshift = True
        _logger.debug(f"check_frameshifts: Detected frameshift due to indel with length_diff={length_diff} (not multiple of 3)")
    else:
        # Otherwise, check if the protein sequences differ in length or content
        is_frameshift = (len(ref_protein_trimmed) != len(edited_protein_trimmed) or 
                        any(r != e for r, e in zip(ref_protein, edited_protein)))
        _logger.debug(f"check_frameshifts: Protein length check - ref: {len(ref_protein_trimmed)}, edited: {len(edited_protein_trimmed)}, is_frameshift: {is_frameshift}")
    
    _logger.debug(f"check_frameshifts: is_indel={is_indel}, length_diff={length_diff}, is_frameshift={is_frameshift}")
    
    # If not a frameshift, check if it's a substitution that introduces a premature stop
    if not is_frameshift and not is_indel and len(ref_cds) == len(edited_cds):
        # This is a substitution, not a frameshift
        _logger.debug("check_frameshifts: Not a frameshift (substitution)")
        return result
    
    # If it's an in-frame indel, it's not a frameshift
    if is_indel and length_diff % 3 == 0 and len(ref_protein) == len(edited_protein):
        _logger.debug("check_frameshifts: In-frame indel, not a frameshift")
        return result
        
    _logger.debug("check_frameshifts: Detected frameshift")
    _logger.debug(f"check_frameshifts: ref_cds: {ref_cds[:30]}...")
    _logger.debug(f"check_frameshifts: edited_cds: {edited_cds[:30]}...")
    _logger.debug(f"check_frameshifts: ref_protein: {ref_protein}")
    _logger.debug(f"check_frameshifts: edited_protein: {edited_protein}")
    
    # This is a frameshift (insertion or deletion not divisible by 3)
    ref_stop_pos = ref_protein.find('*') if '*' in ref_protein else len(ref_protein)
    _logger.debug(f"check_frameshifts: Reference stop position: {ref_stop_pos}")
    
    # Find where the frameshift occurs in the protein sequence
    fs_pos = 0
    for i, (ref_aa, edit_aa) in enumerate(zip(ref_protein, edited_protein)):
        if i >= len(edited_protein) or ref_aa != edit_aa:
            fs_pos = i + 1  # 1-based position
            break
    
    _logger.debug(f"check_frameshifts: Frameshift detected at position {fs_pos}")
    
    # Check for stop gain in the edited protein (both '*' and 'X' are stop codons)
    has_stop = '*' in edited_protein or 'X' in edited_protein
    
    if has_stop:
        # Find the first stop codon (either '*' or 'X')
        stop_pos_star = edited_protein.find('*') if '*' in edited_protein else float('inf')
        stop_pos_X = edited_protein.find('X') if 'X' in edited_protein else float('inf')
        edited_stop_pos = min(stop_pos_star, stop_pos_X)
        
        _logger.debug(f"check_frameshifts: Edited stop at position {edited_stop_pos}, reference stop at {ref_stop_pos}, ref protein length: {len(ref_protein)}")
        
        # For frameshifts, any stop is considered a premature stop
        if is_frameshift:
            is_premature_stop = True
        else:
            # For non-frameshifts, use the standard criteria
            is_premature_stop = (
                (edited_stop_pos < ref_stop_pos) or 
                (ref_stop_pos == len(ref_protein) and has_stop) or
                (edited_stop_pos < len(ref_protein) - 1 and not ref_protein.endswith('*'))
            )
        
        _logger.debug(f"check_frameshifts: is_premature_stop={is_premature_stop} (edited_stop_pos={edited_stop_pos}, ref_stop_pos={ref_stop_pos}, len(ref_protein)={len(ref_protein)})")
        
        if is_premature_stop:
            _logger.debug(f"check_frameshifts: Found stop_gain at position {edited_stop_pos + 1} (frameshift at {fs_pos})")
            
            result.consequence = "stop_gain"
            result.hgvs_p = f"p.{ref_protein[fs_pos-1] if fs_pos > 0 else '?'}{fs_pos}fs*{edited_stop_pos + 1}"
            result.protein_sequence = edited_protein
            return result
    else:
        _logger.debug("check_frameshifts: No stop codon found in edited protein")
    
    # If we get here, it's a regular frameshift without a premature stop
    _logger.debug("check_frameshifts: No premature stop found, returning frameshift_variant")
    result.consequence = "frameshift_variant"
    result.hgvs_p = f"p.{ref_protein[fs_pos-1] if fs_pos > 0 else '?'}{fs_pos}fs"
    result.protein_sequence = edited_protein
    return result

def check_stop_loss(
    ref_cds: str,
    edited_cds: str,
    ref_protein: str,
    edited_protein: str
) -> ConsequenceCheckResult:
    """Check for stop loss variants."""
    result = ConsequenceCheckResult()
    
    # Check for stop loss (must be after stop gain check)
    ref_has_stop = ref_cds[-3:] in STOP_CODONS if len(ref_cds) >= 3 else False
    edited_has_stop = edited_cds[-3:] in STOP_CODONS if len(edited_cds) >= 3 else False
    
    if ref_has_stop and not edited_has_stop:
        # Check if the sequences are the same except for the stop codon
        if ref_cds[:-3] == edited_cds[:-3]:
            result.consequence = "stop_loss"
            ref_aa = translate_sequence(ref_cds[-3:])
            result.hgvs_p = f"p.{ref_aa}{len(edited_protein)}?"
            result.protein_sequence = edited_protein
    
    return result

def check_stop_gain(
    ref_cds: str,
    edited_cds: str,
    ref_protein: str,
    edited_protein: str
) -> ConsequenceCheckResult:
    """Check for stop gains and stop losses from nucleotide changes."""
    result = ConsequenceCheckResult()
    _logger = logging.getLogger(__name__)
    
    # For stop gain/loss, we need to check both sequences
    # First check for stop loss (reference has stop but edited doesn't)
    ref_has_stop = '*' in ref_protein or 'X' in ref_protein
    edited_has_stop = '*' in edited_protein or 'X' in edited_protein
    
    # Check if the reference has a stop codon that's changed to an amino acid in the edited sequence
    if ref_has_stop:
        # Find the position of the stop in the reference
        ref_stop_pos_star = ref_protein.find('*') if '*' in ref_protein else float('inf')
        ref_stop_pos_X = ref_protein.find('X') if 'X' in ref_protein else float('inf')
        ref_stop_pos = min(ref_stop_pos_star, ref_stop_pos_X)
        
        # If we found a stop in the reference
        if ref_stop_pos != float('inf'):
            # Check if the edited sequence doesn't have a stop at the same position
            if (ref_stop_pos >= len(edited_protein) or 
                (edited_protein[ref_stop_pos] != '*' and 
                 edited_protein[ref_stop_pos] != 'X')):
                result.consequence = "stop_loss"
                stop_pos = ref_stop_pos + 1  # 1-based position
                new_aa = edited_protein[ref_stop_pos] if ref_stop_pos < len(edited_protein) else '?'
                
                # If the stop was at the end, use 'ext*' notation
                if ref_stop_pos == len(ref_protein) - 1:
                    result.hgvs_p = f"p.*{stop_pos}ext*"
                else:
                    result.hgvs_p = f"p.*{stop_pos}{new_aa}"
                
                result.protein_sequence = edited_protein
                _logger.debug(f"check_stop_gain: Found stop_lost at position {stop_pos} (stop changed to {new_aa})")
                return result
    
    # Find the first stop in the reference (either '*' or 'X')
    ref_stop_pos_star = ref_protein.find('*') if '*' in ref_protein else float('inf')
    ref_stop_pos_X = ref_protein.find('X') if 'X' in ref_protein else float('inf')
    ref_stop_pos = min(ref_stop_pos_star, ref_stop_pos_X)
    
    # Check for stop loss: reference has a stop that's not in the edited sequence
    # or the stop in the edited sequence is further downstream
    if ref_has_stop:
        if not edited_has_stop:
            # This is a clear stop_lost variant
            if ref_stop_pos != float('inf'):
                result.consequence = "stop_lost"
                # Use the position where the stop was lost (1-based)
                stop_pos = ref_stop_pos + 1
                # Get the new amino acid at that position if it exists
                new_aa = edited_protein[ref_stop_pos] if ref_stop_pos < len(edited_protein) else '?'
                # If the stop was at the end, use 'ext*' notation
                if ref_stop_pos == len(ref_protein) - 1:
                    result.hgvs_p = f"p.*{stop_pos}ext*"
                else:
                    result.hgvs_p = f"p.*{stop_pos}{new_aa}"
                result.protein_sequence = edited_protein
                _logger.debug(f"check_stop_gain: Found stop_lost at position {stop_pos} (stop removed)")
                return result
        else:
            # Both have stops, check if the stop moved downstream
            edited_stop_pos_star = edited_protein.find('*') if '*' in edited_protein else float('inf')
            edited_stop_pos_X = edited_protein.find('X') if 'X' in edited_protein else float('inf')
            edited_stop_pos = min(edited_stop_pos_star, edited_stop_pos_X)
            
            # If the stop moved downstream, it's a stop_lost variant
            if edited_stop_pos > ref_stop_pos and ref_stop_pos != float('inf'):
                result.consequence = "stop_lost"
                stop_pos = ref_stop_pos + 1
                new_aa = edited_protein[ref_stop_pos] if ref_stop_pos < len(edited_protein) else '?'
                result.hgvs_p = f"p.*{stop_pos}{new_aa}"
                result.protein_sequence = edited_protein
                _logger.debug(f"check_stop_gain: Found stop_lost (moved stop) at position {stop_pos}")
                return result
    
    # Also check for the case where the reference has a stop but the edited version has a different stop position
    if ref_has_stop and edited_has_stop:
        ref_stop_pos_star = ref_protein.find('*') if '*' in ref_protein else float('inf')
        ref_stop_pos_X = ref_protein.find('X') if 'X' in ref_protein else float('inf')
        ref_stop_pos = min(ref_stop_pos_star, ref_stop_pos_X)
        
        edited_stop_pos_star = edited_protein.find('*') if '*' in edited_protein else float('inf')
        edited_stop_pos_X = edited_protein.find('X') if 'X' in edited_protein else float('inf')
        edited_stop_pos = min(edited_stop_pos_star, edited_stop_pos_X)
        
        # If the stop moved to a later position, it's a stop_lost variant
        if edited_stop_pos > ref_stop_pos and ref_stop_pos != float('inf'):
            result.consequence = "stop_lost"
            stop_pos = ref_stop_pos + 1
            new_aa = edited_protein[ref_stop_pos] if ref_stop_pos < len(edited_protein) else '?'
            result.hgvs_p = f"p.*{stop_pos}{new_aa}"
            result.protein_sequence = edited_protein
            _logger.debug(f"check_stop_gain: Found stop_lost (moved stop) at position {stop_pos}")
            return result
    
    # Only check for stop gain if lengths are the same (substitution)
    if len(ref_cds) != len(edited_cds):
        _logger.debug("check_stop_gain: Lengths differ, not a simple substitution")
        return result
    
    # Count the number of differences
    diff_positions = [i for i, (ref_nt, edit_nt) in enumerate(zip(ref_cds, edited_cds)) 
                     if ref_nt != edit_nt]
    
    _logger.debug(f"check_stop_gain: Found {len(diff_positions)} differences at positions {diff_positions}")
    _logger.debug(f"check_stop_gain: ref_cds: {ref_cds[:30]}...")
    _logger.debug(f"check_stop_gain: edited_cds: {edited_cds[:30]}...")
    _logger.debug(f"check_stop_gain: ref_protein: {ref_protein}")
    _logger.debug(f"check_stop_gain: edited_protein: {edited_protein}")
    
    # Find the first position where the proteins differ
    first_diff_pos = None
    for i, (ref_aa, edit_aa) in enumerate(zip(ref_protein, edited_protein)):
        if ref_aa != edit_aa:
            first_diff_pos = i
            break
    
    if first_diff_pos is None:
        _logger.debug("check_stop_gain: No differences in protein sequence")
        return result
    
    # Check if the edited protein has a stop codon at the first difference position
    if first_diff_pos is not None and edited_protein[first_diff_pos] == '*':
        # It's a stop gain if the reference doesn't have a stop at this position
        if first_diff_pos >= len(ref_protein) or ref_protein[first_diff_pos] != '*':
            _logger.debug(f"check_stop_gain: Found stop_gain at position {first_diff_pos + 1}")
            
            # Get the reference amino acid at this position
            ref_aa = ref_protein[first_diff_pos] if first_diff_pos < len(ref_protein) else '?'
            
            result.consequence = "stop_gain"
            result.hgvs_p = f"p.{ref_aa}{first_diff_pos + 1}*"
            result.protein_sequence = edited_protein
    # Also check if the edited protein has a stop codon before the reference stop
    elif '*' in edited_protein:
        edited_stop_pos = edited_protein.find('*')
        ref_stop_pos = ref_protein.find('*') if '*' in ref_protein else len(ref_protein)
        
        # It's a stop gain if the stop appears earlier in the edited protein
        # and it's not just the natural stop being moved
        is_premature_stop = (
            (edited_stop_pos < ref_stop_pos) or 
            (ref_stop_pos == -1 and '*' in edited_protein) or
            (edited_stop_pos < len(ref_protein) - 1 and not ref_protein.endswith('*'))
        )
        
        if is_premature_stop:
            # Find the first difference that could have caused the stop
            for pos in range(min(len(ref_protein), len(edited_protein))):
                if ref_protein[pos] != edited_protein[pos]:
                    result.consequence = "stop_gain"
                    result.hgvs_p = f"p.{ref_protein[pos] if pos < len(ref_protein) else '?'}{pos + 1}*"
                    result.protein_sequence = edited_protein
                    break
    
    _logger.debug(f"check_stop_gain: Returning result: {result}")
    return result

def check_missense_synonymous(
    ref_cds: str,
    edited_cds: str,
    ref_protein: str,
    edited_protein: str
) -> ConsequenceCheckResult:
    """Check for missense and synonymous variants."""
    result = ConsequenceCheckResult()
    
    # Only check for missense/synonymous if it's not an indel
    if len(ref_cds) != len(edited_cds):
        return result
    
    # Skip if there's a stop gain (handled by check_stop_gain)
    if '*' in edited_protein:
        return result
    
    # Find the first position where the proteins differ
    for i, (ref_aa, edit_aa) in enumerate(zip(ref_protein, edited_protein)):
        if i >= len(edited_protein):
            break
            
        if ref_aa != edit_aa:
            # Check if it's a stop codon change
            if i == len(ref_protein) - 1 and ref_aa == '*' and edit_aa != '*':
                # This is a stop loss, which we already checked for
                continue
            
            # It's a missense variant
            result.consequence = "missense_variant"
            result.hgvs_p = f"p.{ref_aa}{i+1}{edit_aa}"
            result.protein_sequence = edited_protein
            return result
    
    # If we get here, it's a synonymous variant
    if ref_protein == edited_protein:
        result.consequence = "synonymous_variant"
        result.hgvs_p = "p.="
        result.protein_sequence = edited_protein
    
    return result

def determine_consequence_priority(checks: List[ConsequenceCheckResult]) -> ConsequenceCheckResult:
    """
    Determine the highest priority consequence from a list of check results.
    
    The priority order is:
    1. start_loss
    2. stop_gain
    3. frameshift_variant
    4. stop_loss
    5. in_frame_indel
    6. missense_variant
    7. synonymous_variant
    """
    CONSEQUENCE_PRIORITY = {
        'start_loss': 1,
        'stop_gain': 2,
        'frameshift_variant': 3,
        'stop_loss': 4,
        'in_frame_indel': 5,
        'missense_variant': 6,
        'synonymous_variant': 7
    }
    
    # Filter out negative results
    positive_checks = [c for c in checks if c.is_positive()]
    
    if not positive_checks:
        return ConsequenceCheckResult()
    
    # Return the highest priority consequence
    return min(
        positive_checks,
        key=lambda x: CONSEQUENCE_PRIORITY.get(x.consequence, float('inf'))
    )

def analyze_consequences(
    ref_cdna: str,
    edited_cdna: str,
    config: TranscriptConfig
) -> ProteinOutcome:
    """
    Analyzes the consequences of variants on the protein sequence following standard VEP terminology.
    
    This function performs the following steps:
    1. Extracts CDS from reference and edited sequences
    2. Validates the CDS sequences
    3. Translates both sequences to protein
    4. Runs all consequence checks
    5. Determines the highest priority consequence
    6. Returns a ProteinOutcome with the result
    """
    _logger = logging.getLogger(__name__)
    _logger.setLevel(logging.DEBUG)
    
    # Log input parameters
    _logger.debug("=" * 80)
    _logger.debug("ANALYZE_CONSEQUENCES CALLED")
    _logger.debug("=" * 80)
    _logger.debug(f"Reference cDNA length: {len(ref_cdna) if ref_cdna else 0} bp")
    _logger.debug(f"Edited cDNA length: {len(edited_cdna) if edited_cdna else 0} bp")
    
    # Log first 100bp of each sequence
    _logger.debug(f"Reference cDNA (first 100bp): {ref_cdna[:100] if ref_cdna else 'None'}")
    _logger.debug(f"Edited cDNA (first 100bp): {edited_cdna[:100] if edited_cdna else 'None'}")
    
    # Validate input sequences
    if not ref_cdna or not edited_cdna:
        _logger.error("Empty reference or edited cDNA sequence provided")
        return ProteinOutcome(protein_sequence=None, hgvs_p=None, consequence="no_valid_sequence")
    
    if config.cds_start is None or config.cds_end is None:
        _logger.warning("CDS start or end is None, treating as non-coding transcript")
        return ProteinOutcome(protein_sequence=None, hgvs_p=None, consequence="non_coding_transcript")
    
    # Ensure CDS start and end are within bounds and 1-based
    cds_start = max(1, config.cds_start)
    cds_end = min(config.cds_end, len(ref_cdna))
    
    if cds_start != config.cds_start or cds_end != config.cds_end:
        _logger.warning(f"Adjusted CDS coordinates from [{config.cds_start}-{config.cds_end}] to [{cds_start}-{cds_end}] to fit within sequence length {len(ref_cdna)}")
    
    # Extract CDS sequences (0-based, end-exclusive)
    ref_cds = ref_cdna[cds_start-1:cds_end]
    edited_cds = edited_cdna[cds_start-1:cds_end]
    
    # Log CDS information
    _logger.debug(f"Reference CDS length: {len(ref_cds)}")
    _logger.debug(f"Edited CDS length: {len(edited_cds)}")
    
    # Log CDS lengths for debugging
    if len(ref_cds) % 3 != 0 or len(edited_cds) % 3 != 0:
        _logger.warning(f"CDS lengths not multiples of 3 - ref: {len(ref_cds)}, edited: {len(edited_cds)}")
    
    # For frameshift detection, we need to use the full CDS sequences
    # even if they're not multiples of 3
    is_frameshift = (len(ref_cds) - len(edited_cds)) % 3 != 0
    
    # Create truncated versions for translation if needed
    if not is_frameshift and (len(ref_cds) % 3 != 0 or len(edited_cds) % 3 != 0):
        # For non-frameshift variants, we can truncate to a multiple of 3
        truncated_ref_cds = ref_cds[:-(len(ref_cds) % 3)] if len(ref_cds) % 3 != 0 else ref_cds
        truncated_edited_cds = edited_cds[:-(len(edited_cds) % 3)] if len(edited_cds) % 3 != 0 else edited_cds
        
        if not truncated_ref_cds or not truncated_edited_cds:
            _logger.error("Empty CDS sequence after truncation")
            return ProteinOutcome(protein_sequence=None, hgvs_p=None, consequence="no_valid_CDS")
            
        # Use truncated sequences for translation
        ref_protein = translate_sequence(truncated_ref_cds)
        edited_protein = translate_sequence(truncated_edited_cds)
    else:
        # For frameshifts, we need to use the full sequences
        # even if they're not multiples of 3
        ref_protein = translate_sequence(ref_cds)
        edited_protein = translate_sequence(edited_cds)
    
    _logger.debug(f"Reference protein: {ref_protein}")
    _logger.debug(f"Edited protein:   {edited_protein}")
    
    # First, check for stop_gain/stop_loss as they are high-priority consequences
    stop_gain_result = check_stop_gain(ref_cds, edited_cds, ref_protein, edited_protein)
    _logger.debug(f"analyze_consequences: stop_gain_result = {stop_gain_result}")
    if stop_gain_result.is_positive():
        if stop_gain_result.consequence == 'stop_loss':
            _logger.debug("analyze_consequences: Found stop_loss from check_stop_gain")
            _logger.debug(f"analyze_consequences: stop_gain_result.consequence = {stop_gain_result.consequence}")
            _logger.debug(f"analyze_consequences: stop_gain_result.hgvs_p = {stop_gain_result.hgvs_p}")
            return ProteinOutcome(
                protein_sequence=stop_gain_result.protein_sequence or edited_protein,
                hgvs_p=stop_gain_result.hgvs_p or "p.=?",
                consequence=stop_gain_result.consequence
            )
        elif stop_gain_result.consequence == 'stop_gain':
            _logger.debug("analyze_consequences: Found stop_gain from check_stop_gain")
            return ProteinOutcome(
                protein_sequence=stop_gain_result.protein_sequence or edited_protein,
                hgvs_p=stop_gain_result.hgvs_p or "p.=?",
                consequence=stop_gain_result.consequence
            )
    else:
        _logger.debug("analyze_consequences: No stop_gain or stop_loss found")
    
    # Check for frameshifts next
    frameshift_result = check_frameshifts(ref_cds, edited_cds, ref_protein, edited_protein)
    if frameshift_result.is_positive():
        if frameshift_result.consequence == 'stop_gain':
            _logger.debug("analyze_consequences: Found stop_gain from frameshift")
            return ProteinOutcome(
                protein_sequence=frameshift_result.protein_sequence or edited_protein,
                hgvs_p=frameshift_result.hgvs_p or "p.=?",
                consequence=frameshift_result.consequence
            )
        elif frameshift_result.consequence == 'frameshift_variant':
            _logger.debug("analyze_consequences: Found frameshift_variant")
            return ProteinOutcome(
                protein_sequence=frameshift_result.protein_sequence or edited_protein,
                hgvs_p=frameshift_result.hgvs_p or "p.=?",
                consequence=frameshift_result.consequence
            )
    
    # Run other consequence checks in priority order, using the full CDS sequences
    checks = [
        # Check for start codon loss (highest priority)
        check_start_loss(ref_cds, edited_cds, ref_protein, edited_protein),
        
        # We already checked stop_gain/stop_loss above, so we can skip it here
        
        # Check for in-frame indels
        check_in_frame_indels(ref_cds, edited_cds, ref_protein, edited_protein),
        
        # Check for stop loss
        check_stop_loss(ref_cds, edited_cds, ref_protein, edited_protein),
        
        # Check for missense and synonymous variants (lowest priority)
        check_missense_synonymous(ref_cds, edited_cds, ref_protein, edited_protein)
    ]
    
    # Determine the highest priority consequence
    result = determine_consequence_priority(checks)
    
    # If no consequence was found, return a default outcome
    if not result.is_positive():
        return ProteinOutcome(
            protein_sequence=edited_protein,
            hgvs_p="p.=?",
            consequence="no_sequence_alteration"
        )
    
    # Return the highest priority consequence
    return ProteinOutcome(
        protein_sequence=result.protein_sequence or edited_protein,
        hgvs_p=result.hgvs_p or "p.=?",
        consequence=result.consequence
    )
