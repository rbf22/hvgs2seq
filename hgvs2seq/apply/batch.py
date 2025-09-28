from typing import List, cast
import logging

from ..models import EditPlan
from ..parse import VariantNorm
from .single import apply_single_variant

logger = logging.getLogger(__name__)

def _detect_overlap(v1: VariantNorm, v2: VariantNorm) -> bool:
    """
    Checks if two variants, v1 and v2, have overlapping cDNA coordinates.
    Assumes variants are sorted by start position, so v1.c_start <= v2.c_start.
    """
    return v1.c_end >= v2.c_start

def apply_edit_plan(ref_seq: str, plan: EditPlan) -> List[str]:
    """
    Applies a batch of variants to a reference sequence according to an EditPlan.
    This function handles multiple haplotypes and variant application policies.

    Args:
        ref_seq: The reference cDNA sequence.
        plan: The EditPlan containing variants grouped by haplotype and the policy.

    Returns:
        A list of edited cDNA sequences, one for each haplotype.
    """
    edited_sequences = []

    for i, haplotype_variants in enumerate(plan.haplotypes):
        logger.info(f"Processing haplotype {i+1} with {len(haplotype_variants)} variants.")

        # Sort variants by start position. For deletions, reverse order is safer.
        # For now, we will sort ascending for overlap detection and then reverse for application.
        sorted_variants = sorted(haplotype_variants, key=lambda v: v.c_start)

        if plan.policy == "reject_overlaps":
            for j in range(len(sorted_variants) - 1):
                v1 = sorted_variants[j]
                v2 = sorted_variants[j+1]
                if _detect_overlap(v1, v2):
                    overlap_msg = (
                        f"Overlapping variants detected in haplotype {i+1}: "
                        f"{v1.hgvs_c} (ends at {v1.c_end}) and {v2.hgvs_c} (starts at {v2.c_start}). "
                        "Rejecting this haplotype."
                    )
                    plan.warnings.append(overlap_msg)
                    logger.warning(overlap_msg)
                    # For a strict rejection policy, we might add a placeholder or raise an error.
                    # Here, we'll just return an empty sequence for the failed haplotype.
                    edited_sequences.append("")
                    continue

        # Apply variants in reverse order of start position to handle coordinate shifts correctly.
        variants_to_apply = sorted(haplotype_variants, key=lambda v: v.c_start, reverse=True)

        edited_seq = ref_seq
        for variant in variants_to_apply:
            try:
                edited_seq = apply_single_variant(edited_seq, variant)
            except ValueError as e:
                error_msg = f"Failed to apply variant {variant.hgvs_c} to haplotype {i+1}: {e}"
                plan.warnings.append(error_msg)
                logger.error(error_msg)
                # If a variant fails, we stop processing this haplotype to avoid incorrect results.
                edited_seq = "" # Mark as failed
                break

        edited_sequences.append(edited_seq)

    return edited_sequences