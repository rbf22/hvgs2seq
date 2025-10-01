"""Module for detecting inframe insertions and deletions."""
from typing import Optional
from .base import ConsequenceChecker

class InframeIndelChecker(ConsequenceChecker):
    """Check for inframe insertions and deletions."""

    @classmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for inframe insertions or deletions.

        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence
            edited_protein: Edited protein sequence

        Returns:
            "inframe_insertion" or "inframe_deletion" if detected, None otherwise
        """
        # If either sequence is empty, not an inframe indel
        if not ref_cdna or not edited_cdna:
            return None
            
        # Calculate length difference
        length_diff = len(edited_cdna) - len(ref_cdna)
        
        # Must be a non-zero multiple of 3
        if length_diff == 0 or length_diff % 3 != 0:
            return None
            
        # If we have protein sequences, verify they match the expected pattern for an inframe indel
        if ref_protein and edited_protein:
            # For an inframe insertion, the edited protein should be the same length or longer
            # and should contain the reference protein as a substring (with possible changes)
            if length_diff > 0 and len(edited_protein) < len(ref_protein):
                return None
                
            # For an inframe deletion, the edited protein should be the same length or shorter
            # and should be a substring of the reference protein (with possible changes)
            if length_diff < 0 and len(edited_protein) > len(ref_protein):
                return None
        
        # If we get here and the length difference is a multiple of 3, it's an inframe indel
        if length_diff > 0:
            return "inframe_insertion"
        else:
            return "inframe_deletion"
