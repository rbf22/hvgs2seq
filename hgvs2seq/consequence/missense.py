"""Module for detecting missense variants."""
from typing import Optional
from .base import ConsequenceChecker

class MissenseChecker(ConsequenceChecker):
    """Check for missense variants."""

    @classmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for missense variants.

        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence
            edited_protein: Edited protein sequence

        Returns:
            "missense_variant" if a missense variant is detected, None otherwise
        """
        # If either protein sequence is empty, can't determine missense
        if not ref_protein or not edited_protein:
            return None
            
        # If lengths are different, not a missense (could be frameshift or inframe indel)
        if len(ref_protein) != len(edited_protein):
            return None
            
        # Check for any amino acid differences that aren't stop codons
        for ref_aa, edit_aa in zip(ref_protein, edited_protein):
            # If we find any difference that's not a stop codon, it's a missense
            if ref_aa != edit_aa and ref_aa != '*' and edit_aa != '*':
                return "missense_variant"
                
        return None
