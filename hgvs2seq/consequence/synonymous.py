"""Module for detecting synonymous variants."""
from typing import Optional
from .base import ConsequenceChecker

class SynonymousChecker(ConsequenceChecker):
    """Check for synonymous (silent) variants."""

    @classmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for synonymous variants.

        A variant is considered synonymous if:
        1. The protein sequences are identical
        2. The DNA sequences are different
        3. The variant is in a coding region

        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence
            edited_protein: Edited protein sequence

        Returns:
            "synonymous_variant" if a synonymous variant is detected, None otherwise
        """
        # If either sequence is empty, can't determine if synonymous
        if not ref_cdna or not edited_cdna or not ref_protein or not edited_protein:
            return None
            
        # If sequences are identical, not a variant
        if ref_cdna == edited_cdna:
            return None
            
        # If protein sequences are different, not synonymous
        if ref_protein != edited_protein:
            return None
            
        # If we get here, DNA differs but protein is the same - it's synonymous
        return "synonymous_variant"
