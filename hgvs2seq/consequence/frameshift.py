from typing import Optional
from hgvs2seq.models import TranscriptConfig
from .base import ConsequenceChecker

class FrameshiftChecker(ConsequenceChecker):
    """Check for frameshift variants.
    
    A frameshift is detected when the length difference between the reference
    and edited cDNA is not a multiple of 3, or when the variant causes
    a change in protein length that's not a multiple of 3.
    """

    @classmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for frameshift variants.

        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration (unused in this checker)
            ref_protein: Reference protein sequence (for validation)
            edited_protein: Edited protein sequence (for validation)

        Returns:
            "frameshift_variant" if a frameshift is detected, None otherwise
        """
        # 1. Basic sequence validation - only check ref_cdna, allow empty edited_cdna
        if not ref_cdna:
            return None
            
        # 2. Special case: Empty edited sequence is always a frameshift deletion
        if not edited_cdna:
            return "frameshift_variant"
            
        # 3. Check for frameshift-inducing indels
        # For indels, check if the length change is a multiple of 3
        length_diff = len(edited_cdna) - len(ref_cdna)
        
        # If we have protein sequences, use them to determine if it's a frameshift or inframe indel
        if ref_protein is not None and edited_protein is not None:
            # If protein sequences are the same length but DNA lengths differ, it's an inframe indel
            if len(ref_protein) == len(edited_protein) and length_diff != 0:
                return None
            # If protein sequences differ in length by more than 1, it's a frameshift
            if abs(len(ref_protein) - len(edited_protein)) > 1:
                return "frameshift_variant"
        
        # If we don't have protein sequences, fall back to DNA length check
        # If lengths are the same but we don't have protein sequences, can't determine if it's a frameshift
        if length_diff == 0:
            return None
            
        # For indels, check if the length change is not a multiple of 3
        if length_diff % 3 != 0:
            return "frameshift_variant"
            
        # If we get here, it's an inframe indel (length difference is a multiple of 3)
        # and we either don't have protein sequences or they match our expectations for an inframe indel
        return None
                
        # 5. If we get here, it's not a frameshift
        return None