"""Module for detecting start loss variants."""
from typing import Optional
from .base import ConsequenceChecker

class StartLossChecker(ConsequenceChecker):
    """Check for start loss variants."""

    @classmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for start loss variants.

        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence
            edited_protein: Edited protein sequence

        Returns:
            "start_lost" if a start loss variant is detected, None otherwise
        """
        # If either protein sequence is empty, can't determine start loss
        if not ref_protein or not edited_protein:
            return None
            
        # If reference doesn't start with 'M', can't have start loss
        if not ref_protein.startswith('M'):
            return None
            
        # Check if the start codon is lost in the edited protein
        if not edited_protein.startswith('M'):
            return "start_lost"
            
        # Check if the start codon was changed to a different amino acid
        if ref_protein.startswith('M') and not edited_protein.startswith('M'):
            return "start_lost"
            
        return None
