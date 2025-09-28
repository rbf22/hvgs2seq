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
    start_idx = variant.c_start - 1
    end_idx = variant.c_end

    if start_idx < 0 or end_idx > len(ref_seq):
        raise ValueError(
            f"Variant coordinates ({variant.c_start}, {variant.c_end}) are outside "
            f"the bounds of the reference sequence (length {len(ref_seq)})."
        )

    prefix = ref_seq[:start_idx]

    if variant.kind == 'ins':
        insertion_idx = variant.c_start
        prefix = ref_seq[:insertion_idx]
        suffix = ref_seq[insertion_idx:]
        return prefix + variant.alt + suffix

    suffix = ref_seq[end_idx:]
    region_to_change = ref_seq[start_idx:end_idx]

    if variant.kind in ("sub", "mnv", "delins"):
        return prefix + variant.alt + suffix

    elif variant.kind == "del":
        return prefix + suffix

    elif variant.kind == "dup":
        # The result is prefix + original_region + duplicated_region + suffix
        return prefix + region_to_change + region_to_change + suffix

    elif variant.kind == "inv":
        # Inversion: replace the affected region with its reverse complement
        inverted_part = _reverse_complement(region_to_change)
        return prefix + inverted_part + suffix

    else:
        _logger.error(f"Unsupported variant kind: {variant.kind}")
        raise ValueError(f"Applying variant of kind '{variant.kind}' is not supported.")