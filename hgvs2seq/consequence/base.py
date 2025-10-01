from typing import Optional, Dict, Any, List
from abc import ABC, abstractmethod
class ConsequenceChecker(ABC):
    """Abstract base class for all consequence checkers."""
    
    @staticmethod
    def translate_sequence(sequence: str) -> str:
        """Translate a DNA sequence to a protein sequence, starting from the first start codon (ATG).
        
        Args:
            sequence: DNA sequence to translate
            
        Returns:
            Translated protein sequence, or empty string if no start codon found
        """
        # Find the first start codon
        start_pos = sequence.upper().find('ATG')
        if start_pos == -1:
            return ""  # No start codon found
            
        # Only translate from the start codon onwards
        coding_sequence = sequence[start_pos:]
        
        codon_table = {
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
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
        
        protein = []
        for i in range(0, len(coding_sequence) - 2, 3):
            codon = coding_sequence[i:i+3].upper()
            if len(codon) == 3:
                aa = codon_table.get(codon, 'X')
                if aa == '*':  # Stop codon
                    break
                protein.append(aa)
        return ''.join(protein)

    @classmethod
    @abstractmethod
    def check(
        cls, 
        ref_cdna: str, 
        edited_cdna: str, 
        config: 'TranscriptConfig',
        ref_protein: str,
        edited_protein: str
    ) -> Optional[str]:
        """Check for a specific consequence.
        
        Args:
            ref_cdna: Reference cDNA sequence
            edited_cdna: Edited cDNA sequence
            config: Transcript configuration
            ref_protein: Reference protein sequence
            edited_protein: Edited protein sequence
            
        Returns:
            Consequence string if the consequence is found, None otherwise
        """
        pass

# Define the severity order of consequences (from most severe to least)
CONSEQUENCE_SEVERITY = [
    'start_lost',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'inframe_deletion',
    'inframe_insertion',
    'missense_variant',
    'synonymous_variant'
]

def get_most_severe_consequence(results: List[str]) -> Optional[str]:
    """Return the most severe consequence from a list of results.
    
    Args:
        results: List of consequence strings
        
    Returns:
        The most severe consequence, or None if the input list is empty
    """
    if not results:
        return None
        
    # Find the consequence with the highest severity (lowest index in CONSEQUENCE_SEVERITY)
    min_index = len(CONSEQUENCE_SEVERITY)
    most_severe = None
    
    for consequence in results:
        try:
            idx = CONSEQUENCE_SEVERITY.index(consequence)
            if idx < min_index:
                min_index = idx
                most_severe = consequence
        except ValueError:
            # If the consequence is not in our severity list, keep the first one we found
            if most_severe is None:
                most_severe = consequence
                
    return most_severe
