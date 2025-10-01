"""Module for detecting stop gain variants."""
from typing import Optional
from .base import ConsequenceChecker

class StopGainChecker(ConsequenceChecker):
    """Check for stop gain (nonsense) variants."""

    @classmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for stop gain variants.

        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence
            edited_protein: Edited protein sequence

        Returns:
            "stop_gained" if a stop gain variant is detected, None otherwise
        """
        # If either protein sequence is empty, can't determine stop gain
        if not ref_protein or not edited_protein:
            return None
            
        # If lengths are different, not a simple stop gain (could be frameshift)
        if len(ref_protein) != len(edited_protein):
            return None
            
        # Check if the edited protein has a stop codon and the reference doesn't have one at the same position
        if '*' in edited_protein:
            # Get the position of the first stop in the edited protein
            edited_stop_pos = edited_protein.index('*')
            
            # If the reference doesn't have a stop at this position, it's a stop gain
            if edited_stop_pos >= len(ref_protein) or ref_protein[edited_stop_pos] != '*':
                return "stop_gained"
        
        # Special case: If the edited protein ends with a stop but is shorter than the reference
        # and the reference doesn't have a stop at that position
        if edited_protein.endswith('*') and len(edited_protein) <= len(ref_protein):
            if not ref_protein.endswith('*') or len(edited_protein) < len(ref_protein.rstrip('*')):
                return "stop_gained"
                
        return None
