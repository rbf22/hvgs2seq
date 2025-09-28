from typing import List, TYPE_CHECKING
from ..parse import VariantNorm
from .single import apply_single_variant

def apply_batch(ref_seq: str, variants: List[VariantNorm]) -> str:
    """
    Applies a batch of normalized variants to a reference cDNA sequence.

    To prevent issues with coordinates shifting after an indel, variants are
    applied in reverse order of their starting position.

    Args:
        ref_seq: The reference cDNA sequence.
        variants: A list of VariantNorm objects to apply.

    Returns:
        The final modified cDNA sequence after all variants have been applied.
    """
    # Sort variants by start position in descending order.
    # This is crucial to ensure that applying one variant doesn't invalidate the
    # coordinates of subsequent variants. For example, a deletion at the start
    # of the sequence would shift all other coordinates. By working from the
    # end of the sequence to the beginning, we avoid this problem.
    sorted_variants = sorted(variants, key=lambda v: v.c_start, reverse=True)

    edited_seq = ref_seq
    for variant in sorted_variants:
        edited_seq = apply_single_variant(edited_seq, variant)

    return edited_seq