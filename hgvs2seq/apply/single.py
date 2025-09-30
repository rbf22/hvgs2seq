"""
Applies a single normalized variant to a reference sequence.
This is the core engine for sequence manipulation.
"""
import logging
from ..models import VariantNorm

_logger = logging.getLogger(__name__)

def _reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATCGN", "TAGCN")
    return seq.upper().translate(complement_map)[::-1]

def apply_single_variant(ref_seq: str, variant: VariantNorm) -> str:
    """
    Applies a single normalized variant to a reference cDNA sequence.

    Args:
        ref_seq: The reference cDNA sequence.
        variant: A VariantNorm object describing the change.

    Returns:
        The edited cDNA sequence.

    Raises:
        ValueError: If the variant type is unsupported or coordinates are invalid.
    """
    _logger.debug(f"Applying variant {variant.hgvs_c} to sequence of length {len(ref_seq)}")

    # Convert from 1-based HGVS coordinates to 0-based Python indices
    # For inversions, we need to adjust the end index to be exclusive
    start_idx = variant.start - 1
    end_idx = variant.end

    if start_idx < 0 or end_idx > len(ref_seq):
        raise ValueError(
            f"Variant coordinates ({variant.start}, {variant.end}) are outside "
            f"the bounds of the reference sequence (length {len(ref_seq)})."
        )

    prefix = ref_seq[:start_idx]

    if variant.variant_type == 'ins':
        insertion_idx = variant.start
        prefix = ref_seq[:insertion_idx]
        suffix = ref_seq[insertion_idx:]
        return prefix + variant.alt + suffix

    suffix = ref_seq[end_idx:]
    region_to_change = ref_seq[start_idx:end_idx]

    if variant.variant_type in ("sub", "mnv", "delins"):
        return prefix + variant.alt + suffix

    elif variant.variant_type == "del":
        return prefix + suffix

    elif variant.variant_type == "dup":
        # The result is prefix + original_region + duplicated_region + suffix
        return prefix + region_to_change + region_to_change + suffix

    elif variant.variant_type == "inv":
        # Inversion: replace the affected region with its reverse complement
        # For c.2_4inv, we want to invert positions 2-4 (1-based), which is indices 1-4 (0-based, end-exclusive)
        # The region_to_change is ref_seq[1:4] for c.2_4inv ("TGC")
        # Reverse complement of "TGC" is "GCA"
        # So the result should be ref_seq[0:1] + "GCA" + ref_seq[4:]
        # Which is 'A' + 'GCA' + 'ACGTGCAT' = 'AGCAACGTGCAT'
        
        # Get the region to invert (1-based to 0-based, end-exclusive)
        # For c.2_4inv, we want to invert positions 2-4 (1-based)
        start = variant.start - 1  # 2 -> 1 (0-based)
        end = variant.end          # 4 (end-exclusive)
        
        # Get the region to invert and reverse complement it
        region_to_invert = ref_seq[start:end]
        inverted_region = _reverse_complement(region_to_invert)
        
        # Debug prints
        print(f"Reference sequence: {ref_seq}")
        print(f"Inverting region: {region_to_invert} (indices {start}:{end})")
        print(f"Reverse complement: {inverted_region}")
        print(f"Prefix: {ref_seq[:start]}")
        print(f"Suffix: {ref_seq[end:]}")
        
        # Combine and return
        # The prefix is everything before the inverted region
        # The suffix is everything after the inverted region
        # For c.2_4inv, this should be:
        # - ref_seq[:1] = 'A' (first base)
        # - inverted_region = 'GCA' (reverse complement of 'TGC')
        # - ref_seq[4:] = 'ACGTGCAT' (everything after position 4, which is index 4)
        # Result: 'A' + 'GCA' + 'ACGTGCAT' = 'AGCAACGTGCAT'
        result = ref_seq[:start] + inverted_region + ref_seq[end:]
        print(f"Result: {result}")
        print(f"Expected: AGCAACGTGCAT")
        
        # For debugging, let's see the actual positions being used
        print(f"Start position: {start}, End position: {end}")
        print(f"Region to invert: {ref_seq[start:end]}")
        print(f"Result length: {len(result)}, Expected length: {len('AGCAACGTGCAT')}")
        
        return result

    else:
        _logger.error(f"Unsupported variant kind: {variant.kind}")
        raise ValueError(f"Applying variant of kind '{variant.kind}' is not supported.")