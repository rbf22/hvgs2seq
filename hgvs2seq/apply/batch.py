"""
Applies a batch of variants to a reference sequence, handling phasing and ordering.
"""
import logging
from typing import List, Dict
from collections import defaultdict
from ..models import VariantNorm
from .single import apply_single_variant

_logger = logging.getLogger(__name__)

def apply_variants_in_batch(
    ref_seq: str,
    variants: List[VariantNorm],
    policy: str = "order_by_pos"
) -> Dict[int, str]:
    """
    Applies a batch of normalized variants to a reference sequence,
    handling phasing into separate haplotypes.

    Args:
        ref_seq: The reference cDNA sequence.
        variants: A list of VariantNorm objects to apply.
        policy: The policy for handling multiple variants. Currently, only
                "order_by_pos" (applying in reverse order) is implemented.

    Returns:
        A dictionary mapping a haplotype ID to its edited sequence.
        Unphased variants are grouped into haplotype 0.
    """
    if policy != "order_by_pos":
        _logger.warning(f"Policy '{policy}' is not fully supported. "
                        "Defaulting to 'order_by_pos' behavior.")

    # 1. Group variants by phase_group (haplotype)
    haplotypes = defaultdict(list)
    for var in variants:
        # Treat unphased variants (phase_group=None) as a single haplotype (group 0)
        phase_group = var.phase_group if var.phase_group is not None else 0
        haplotypes[phase_group].append(var)

    _logger.info(f"Applying variants across {len(haplotypes)} haplotype(s).")

    edited_sequences = {}

    # 2. For each haplotype, apply its variants to a fresh copy of the ref_seq
    for phase_id, hap_variants in haplotypes.items():
        _logger.info(f"Processing haplotype {phase_id} with {len(hap_variants)} variants.")

        current_seq = ref_seq

        # 3. Sort variants in reverse order of position. This is critical to
        # prevent earlier edits from invalidating the coordinates of later edits.
        sorted_variants = sorted(
            hap_variants,
            key=lambda v: (v.start, v.end),
            reverse=True
        )

        # 4. Apply the sorted variants sequentially
        for variant in sorted_variants:
            try:
                current_seq = apply_single_variant(current_seq, variant)
            except ValueError as e:
                _logger.error(f"Failed to apply variant {variant.hgvs_c} "
                              f"to haplotype {phase_id}: {e}")
                # Re-raise to halt execution on an invalid variant application
                raise

        edited_sequences[phase_id] = current_seq
        _logger.info(f"Finished processing haplotype {phase_id}.")

    return edited_sequences