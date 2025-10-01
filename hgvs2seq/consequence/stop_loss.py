"""Module for detecting stop loss variants."""
from typing import Optional
from .base import ConsequenceChecker

class StopLossChecker(ConsequenceChecker):
    """Check for stop loss variants."""

    @classmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for stop loss variants.

        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence
            edited_protein: Edited protein sequence

        Returns:
            "stop_lost" if a stop loss variant is detected, None otherwise
        """
        # If either protein sequence is empty, can't determine stop loss
        if not ref_protein or not edited_protein:
            return None
            
        # If reference doesn't end with a stop codon, can't have stop loss
        if not ref_protein.endswith('*'):
            return None
            
        # Check if the stop codon is lost in the edited protein
        if not edited_protein.endswith('*'):
            return "stop_lost"
        
        # Also check if the stop codon is changed to an amino acid
        # This handles cases where the stop is still present but at a different position
        elif ref_protein.endswith('*') and edited_protein.endswith('*'):
            # If the stop position has changed, it's a stop loss
            if len(edited_protein) > len(ref_protein):
                return "stop_lost"
        
        # If we have a stop in the reference but not in the edited sequence
        elif ref_protein.endswith('*') and not edited_protein.endswith('*'):
            return "stop_lost"
            
        return None
