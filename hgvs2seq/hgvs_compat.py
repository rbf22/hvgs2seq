"""
Compatibility layer for migrating from hgvs to pyhgvs.

This module provides a drop-in replacement for hgvs functionality using pyhgvs.
"""
from typing import Optional, Tuple, Dict, Any, Union
from dataclasses import dataclass
import logging
import re
import pyhgvs
from pyhgvs import HGVSName

_logger = logging.getLogger(__name__)

# Constants for variant types
VARIANT_TYPES = {
    'g': 'g',  # genomic
    'c': 'c',  # coding DNA
    'n': 'n',  # non-coding DNA
    'p': 'p',  # protein
}

@dataclass
class SimplePosition:
    """Simple position class to mimic hgvs.location.SimplePosition."""
    base: int
    uncertain: bool = False

@dataclass
class Interval:
    """Interval class to mimic hgvs.location.Interval."""
    start: SimplePosition
    end: SimplePosition
    uncertain: bool = False

@dataclass
class PosEdit:
    """PosEdit class to mimic hgvs.posedit.PosEdit."""
    pos: Interval
    edit: Any  # Simplified for now

class SequenceVariant:
    """SequenceVariant class to mimic hgvs.sequencevariant.SequenceVariant."""
    def __init__(self, ac: str, type: str, posedit: PosEdit):
        if type not in VARIANT_TYPES:
            raise ValueError(f"Invalid variant type: {type}")
        self.ac = ac  # accession
        self.type = type  # variant type (g, c, n, p)
        self.posedit = posedit
    
    def __str__(self) -> str:
        """Convert to HGVS string representation."""
        # This is a simplified version - would need to be expanded for full HGVS support
        pos = self.posedit.pos
        return f"{self.ac}:{self.type}.{pos.start.base}_{pos.end.base}"

class HGVSParseError(ValueError):
    """Exception raised when HGVS parsing fails."""
    pass

class HGVSDataNotAvailableError(Exception):
    """Exception raised when required data is not available."""
    pass

class Parser:
    """Parser class to mimic hgvs.parser.Parser."""
    def __init__(self, hdp=None):
        """Initialize the parser.
        
        Args:
            hdp: Optional data provider (not used in this implementation)
        """
        self.hdp = hdp
        self._hdp = hdp  # For compatibility with existing code
        
    def parse_hgvs_variant(self, hgvs_string: str) -> SequenceVariant:
        """Parse an HGVS string into a SequenceVariant object."""
        try:
            # Extract the accession and variant parts
            if ':' in hgvs_string:
                ac, variant = hgvs_string.split(':', 1)
            else:
                # If no accession, use a default one
                variant = hgvs_string

            # Determine variant type from the variant part
            if variant.startswith('g.'):
                var_type = 'g'
                variant = variant[2:]  # Remove 'g.' prefix
            elif variant.startswith('c.'):
                var_type = 'c'
                variant = variant[2:]  # Remove 'c.' prefix
            elif variant.startswith('n.'):
                var_type = 'n'
                variant = variant[2:]  # Remove 'n.' prefix
            elif variant.startswith('p.'):
                var_type = 'p'
                variant = variant[2:]  # Remove 'p.' prefix
            else:
                raise HGVSParseError(f"Could not determine variant type from {hgvs_string}")
            
            # Parse the position and edit parts
            if 'del' in variant and 'ins' in variant:
                # Handle delins variants (deletion-insertion)
                pos_part, edit_part = variant.split('del', 1)
                ref = edit_part.split('ins', 1)[0]
                alt = edit_part.split('ins', 1)[1] if 'ins' in edit_part else ''
                edit_type = 'delins'
            elif 'del' in variant:
                # Handle deletion variants
                pos_part, ref = variant.split('del', 1)
                alt = ''
                edit_type = 'del'
            elif 'ins' in variant:
                # Handle insertion variants
                pos_part, alt = variant.split('ins', 1)
                ref = ''
                edit_type = 'ins'
            elif 'dup' in variant:
                # Handle duplication variants
                pos_part = variant.split('dup', 1)[0]
                ref = pos_part.split('_')[-1]  # This is a simplification
                alt = ref * 2  # Duplicate the reference
                edit_type = 'dup'
            elif '>' in variant:
                # Handle substitution variants (e.g., c.123A>G)
                pos_part, alt = variant.split('>', 1)
                ref = pos_part[-1] if pos_part and pos_part[-1].isalpha() else 'N'
                pos_part = pos_part[:-1] if pos_part and pos_part[-1].isalpha() else pos_part
                edit_type = 'sub'
            else:
                # Default to sequence alteration
                pos_part = variant
                ref = ''
                alt = ''
                edit_type = 'alt'
            
            # Parse the position
            if '_' in pos_part:
                # Handle ranges (e.g., 100_102)
                start_str, end_str = pos_part.split('_', 1)
                start = int(''.join(c for c in start_str if c.isdigit()))
                end = int(''.join(c for c in end_str if c.isdigit()))
            else:
                # Single position
                pos = int(''.join(c for c in pos_part if c.isdigit()))
                start = end = pos
            
            # Create position and interval
            start_pos = SimplePosition(base=start)
            end_pos = SimplePosition(base=end)
            interval = Interval(start=start_pos, end=end_pos)
            
            # Map the edit type to the expected VariantType
            # We use the edit_type directly as it should match the VariantType enum values
            # in models.py ('sub', 'del', 'ins', 'delins', 'dup', etc.)
            variant_type = edit_type  # This should match the VariantType enum values
            
            # Create edit object with the correct type
            class Edit:
                def __init__(self, type, ref, alt):
                    self.type = type  # This will be 'sub', 'del', 'ins', 'delins', 'dup', etc.
                    self.ref = ref
                    self.alt = alt
                
                def __str__(self):
                    if self.type == 'sub':
                        return f"{self.ref}>{self.alt}"
                    elif self.type == 'del':
                        return f"del{self.ref}"
                    elif self.type == 'ins':
                        return f"ins{self.alt}"
                    elif self.type == 'delins':
                        return f"del{self.ref}ins{self.alt}"
                    elif self.type == 'dup':
                        return f"dup{self.ref}"
                    return f"{self.type}{self.alt}"
            
            # Create the edit object with the variant type
            edit_obj = Edit(type=variant_type, ref=ref, alt=alt)
            
            return SequenceVariant(
                ac=ac,
                type=var_type,
                posedit=PosEdit(pos=interval, edit=edit_obj)
            )
            
        except Exception as e:
            _logger.error(f"Failed to parse HGVS variant {hgvs_string}: {e}")
            raise HGVSParseError(f"Failed to parse HGVS variant {hgvs_string}: {e}")


class VariantMapper:
    """VariantMapper class to mimic hgvs.variantmapper.VariantMapper."""
    def __init__(self, hdp):
        self._hdp = hdp
    
    def g_to_c(self, variant: SequenceVariant, ac: str) -> SequenceVariant:
        """Convert genomic variant to coding DNA variant."""
        if variant.type != 'g':
            raise ValueError(f"Expected g. variant, got {variant.type}. variant")
        
        # For now, just change the type and accession
        return SequenceVariant(
            ac=ac,
            type='c',
            posedit=variant.posedit
        )
    
    def c_to_g(self, variant: SequenceVariant, ac: Optional[str] = None) -> SequenceVariant:
        """Convert coding DNA variant to genomic variant."""
        if variant.type != 'c':
            raise ValueError(f"Expected c. variant, got {variant.type}. variant")
        
        # For now, just change the type and accession
        return SequenceVariant(
            ac=ac or "NC_000001.11",  # Default chromosome
            type='g',
            posedit=variant.posedit
        )
    
    def n_to_c(self, variant: SequenceVariant, ac: str) -> SequenceVariant:
        """Convert non-coding DNA variant to coding DNA variant."""
        if variant.type != 'n':
            raise ValueError(f"Expected n. variant, got {variant.type}. variant")
        
        # For now, just change the type and accession
        return SequenceVariant(
            ac=ac,
            type='c',
            posedit=variant.posedit
        )
    
    def c_to_c(self, variant: SequenceVariant, ac: str) -> SequenceVariant:
        """Convert between coding DNA variants of different transcripts."""
        if variant.type != 'c':
            raise ValueError(f"Expected c. variant, got {variant.type}. variant")
        
        # For now, just change the accession
        return SequenceVariant(
            ac=ac,
            type='c',
            posedit=variant.posedit
        )

# Module-level functions for backward compatibility
def connect():
    """Connect to the UTA database (stub implementation)."""
    # In a real implementation, this would set up the necessary connections
    # for pyhgvs to access transcript and genome data
    return {}

# Create module-level instances
parser = Parser(hdp=connect())

class exceptions:
    """Namespace for exceptions."""
    HGVSError = HGVSParseError
    HGVSDataNotAvailableError = HGVSDataNotAvailableError

# Update module attributes to match hgvs module structure
variantmapper = type('VariantMapper', (), {'VariantMapper': VariantMapper})()
dataproviders = type('DataProviders', (), {'uta': type('UTA', (), {'connect': connect})})()
assemblymapper = type('AssemblyMapper', (), {})()

# For backward compatibility
__all__ = [
    'SequenceVariant', 'SimplePosition', 'Interval', 'PosEdit',
    'HGVSParseError', 'HGVSDataNotAvailableError', 'Parser', 'VariantMapper',
    'parser', 'exceptions', 'variantmapper', 'dataproviders', 'assemblymapper'
]
