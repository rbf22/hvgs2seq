from typing import List, Optional, Type
from .base import ConsequenceChecker, get_most_severe_consequence
from .start_loss import StartLossChecker
from .stop_gain import StopGainChecker
from .stop_loss import StopLossChecker
from .inframe_indel import InframeIndelChecker
from .frameshift import FrameshiftChecker
from .synonymous import SynonymousChecker
from .missense import MissenseChecker

class ConsequenceAnalyzer:
    """Analyzes variants to determine their most severe consequence."""
    
    def __init__(self, checkers: Optional[List[Type[ConsequenceChecker]]] = None):
        """Initialize with a list of consequence checkers.
        
        Args:
            checkers: List of consequence checker classes to use. If None, uses all available checkers.
        """
        if checkers is None:
            self.checkers = [
                StartLossChecker,
                StopGainChecker,
                StopLossChecker,
                # Check for inframe indels before frameshifts
                InframeIndelChecker,
                FrameshiftChecker,
                SynonymousChecker,
                MissenseChecker
            ]
        else:
            self.checkers = checkers
    
    def analyze(
        self,
        ref_cdna: str,
        edited_cdna: str,
        config: 'TranscriptConfig',
        ref_protein: str = "",
        edited_protein: str = ""
    ) -> Optional[str]:
        """Analyze the variant and return the most severe consequence.
        
        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence (optional)
            edited_protein: Edited protein sequence (optional)
            
        Returns:
            The most severe consequence as a string, or None if no consequences found
        """
        results: List[str] = []
        
        # Calculate length difference and check if it's a frameshift or inframe indel
        length_diff = len(edited_cdna) - len(ref_cdna)
        is_frameshift = length_diff % 3 != 0
        is_inframe_indel = length_diff != 0 and not is_frameshift
        
        # If we have protein sequences, use them to determine if it's a frameshift or inframe indel
        if ref_protein and edited_protein:
            # If protein sequences are the same length but DNA lengths differ, it's an inframe indel
            if len(ref_protein) == len(edited_protein) and length_diff != 0:
                is_frameshift = False
                is_inframe_indel = True
            # If protein sequences differ in length by more than 1, it's a frameshift
            elif abs(len(ref_protein) - len(edited_protein)) > 1:
                is_frameshift = True
                is_inframe_indel = False
        
        for checker in self.checkers:
            try:
                # Skip inframe indel checker if it's a frameshift
                if is_frameshift and checker.__name__ == "InframeIndelChecker":
                    continue
                    
                # Skip frameshift checker if it's an inframe indel or no length change
                if (is_inframe_indel or length_diff == 0) and checker.__name__ == "FrameshiftChecker":
                    continue
                    
                # Skip synonymous checker if the protein sequences are different or empty
                if (not ref_protein or not edited_protein or ref_protein != edited_protein) and checker.__name__ == "SynonymousChecker":
                    continue
                    
                result = checker.check(
                    ref_cdna=ref_cdna,
                    edited_cdna=edited_cdna,
                    config=config,
                    ref_protein=ref_protein,
                    edited_protein=edited_protein
                )
                if result is not None:
                    results.append(result)
            except Exception as e:
                print(f"Error in {checker.__name__}: {e}")
        
        return get_most_severe_consequence(results)
