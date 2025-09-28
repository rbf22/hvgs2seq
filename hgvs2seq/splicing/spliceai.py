"""
Handles annotation of variants that may impact splicing.
Initially, this module will not run the SpliceAI tool itself, but will
implement the logic to flag variants in or near canonical splice sites.
"""
import logging
from typing import List, Optional
from ..models import VariantNorm, TranscriptConfig

_logger = logging.getLogger(__name__)

# Define the windows for splice site annotation
CANONICAL_DONOR_WINDOW = {1, 2}
CANONICAL_ACCEPTOR_WINDOW = {-1, -2}

def annotate_splicing_impact(
    variant: VariantNorm,
    config: TranscriptConfig
) -> Optional[str]:
    """
    Analyzes a single variant to see if it falls within a canonical splice region.

    Args:
        variant: The VariantNorm object to analyze.
        config: The TranscriptConfig containing exon boundary information.

    Returns:
        A string describing the potential splice impact, or None if not applicable.
    """
    # Splicing analysis is only relevant for variants with intronic offsets
    if variant.c_start_offset == 0 and variant.c_end_offset == 0:
        return None

    # Calculate the cDNA coordinates of exon boundaries
    exon_lengths = [(end - start + 1) for start, end in config.exons]
    exon_end_positions_c = {sum(exon_lengths[:i+1]) for i in range(len(exon_lengths))}
    # The start of the first exon is 1. Subsequent starts are after the previous end.
    exon_start_positions_c = {1}
    cumulative_length = 0
    for length in exon_lengths[:-1]:
        cumulative_length += length
        exon_start_positions_c.add(cumulative_length + 1)

    # Check for donor site variants (e.g., c.123+1G>A)
    # These occur after an exon ends.
    if variant.c_start in exon_end_positions_c and variant.c_start_offset > 0:
        if variant.c_start_offset in CANONICAL_DONOR_WINDOW:
            return f"canonical_donor_site_variant (at c.{variant.c_start}+{variant.c_start_offset})"
        else:
            return f"donor_region_variant (at c.{variant.c_start}+{variant.c_start_offset})"

    # Check for acceptor site variants (e.g., c.124-2A>G)
    # These occur before an exon starts.
    if variant.c_start in exon_start_positions_c and variant.c_start_offset < 0:
        if variant.c_start_offset in CANONICAL_ACCEPTOR_WINDOW:
            return f"canonical_acceptor_site_variant (at c.{variant.c_start}{variant.c_start_offset})"
        else:
            return f"acceptor_region_variant (at c.{variant.c_start}{variant.c_start_offset})"

    return "intronic_variant"


def analyze_all_variants_for_splicing(
    variants: List[VariantNorm],
    config: TranscriptConfig
) -> List[str]:
    """
    Analyzes a list of variants and returns annotations for those with
    potential splicing impact.

    Args:
        variants: A list of VariantNorm objects.
        config: The TranscriptConfig.

    Returns:
        A list of string annotations for relevant variants.
    """
    annotations = []
    for var in variants:
        annotation = annotate_splicing_impact(var, config)
        if annotation:
            annotations.append(f"Variant {var.hgvs_c}: {annotation}")

    _logger.info(f"Found {len(annotations)} variants with potential splicing impact.")
    return annotations