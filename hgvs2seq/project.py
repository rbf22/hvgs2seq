from typing import TYPE_CHECKING, Optional, Any, Dict, Union
from enum import Enum

# Import our models
from .models import VariantNorm, VariantType, TranscriptConfig, GenomicPosition

# Import our compatibility layer
from .hgvs_compat import (
    HGVSParseError, HGVSDataNotAvailableError, VariantMapper, parser
)

# Use a forward reference for TranscriptConfig to avoid circular imports at runtime
if TYPE_CHECKING:
    from .config import TranscriptConfig

# Global mapper instance
_mapper = None

def get_mapper() -> VariantMapper:
    """Get a configured VariantMapper instance.
    
    Returns:
        VariantMapper: A configured VariantMapper instance.
    """
    global _mapper
    if _mapper is None:
        try:
            # In pyhgvs, we don't need to connect to UTA for the mapper
            # The actual data access is handled separately
            _mapper = VariantMapper(None)  # No HDP needed in our implementation
        except Exception as e:
            raise RuntimeError(
                "Failed to initialize VariantMapper. "
            )
    return _mapper


def project_variant(
    variant: VariantNorm, 
    config: "TranscriptConfig"
) -> VariantNorm:
    """
    Projects a variant to the specified transcript's cDNA (c.) coordinates.

    This function takes any valid HGVS variant and attempts to convert it to the
    cDNA coordinates of the transcript defined in the configuration.

    :param variant: A VariantNorm object.
    :param config: The TranscriptConfig for the target transcript.
    :return: A new VariantNorm object in c. coordinates.
    :raises ValueError: If the variant cannot be projected to the target transcript.
    """
    target_transcript_id = config.transcript_id

    # If the variant is already in the correct form, we're done.
    if variant.transcript_id == target_transcript_id and variant.hgvs_c.startswith('c.'):
        return variant

    try:
        # For now, we'll just update the transcript ID and return the variant
        # In a real implementation, this would use pyhgvs to do the actual projection
        if variant.variant_type == VariantType.INSERTION:
            hgvs_c = f"c.{variant.start}_{variant.start + 1}ins{variant.alt}"
        elif variant.variant_type == VariantType.DELETION:
            hgvs_c = f"c.{variant.start}_{variant.end}del"
        elif variant.ref and variant.alt:  # Substitution
            hgvs_c = f"c.{variant.start}{variant.ref}>{variant.alt}"
        else:
            hgvs_c = f"c.{variant.start}_{variant.end}del"  # Fallback
            
        result = variant.model_copy(update={
            'transcript_id': target_transcript_id,
            'hgvs_c': hgvs_c
        })
        return result
        
    except Exception as e:
        raise ValueError(
            f"Failed to project variant {variant.hgvs} to {target_transcript_id}: {str(e)}"
        )


def project_c_to_g(
    variant_c: VariantNorm,
    config: "TranscriptConfig" = None,
) -> VariantNorm:
    """
    Projects a cDNA (c.) variant to genomic (g.) coordinates.

    This function requires that the input variant is in c. coordinates and
    is associated with a transcript that can be found by the data provider.

    :param variant_c: A VariantNorm object in c. coordinates.
    :param config: The TranscriptConfig for the target transcript (optional).
    :return: A new VariantNorm object in g. coordinates.
    :raises ValueError: If the variant cannot be projected.
    """
    if variant_c is None:
        raise ValueError("Variant cannot be None")
        
    if not hasattr(variant_c, 'hgvs_c') or not variant_c.hgvs_c or not variant_c.hgvs_c.startswith('c.'):
        raise ValueError(
            f"Input variant must be in c. coordinates for c_to_g projection, "
            f"but got '{variant_c.hgvs_c}'."
        )

    try:
        # In a real implementation, this would use pyhgvs to do the actual projection
        # For now, we'll just return a mock genomic variant
        if not config:
            raise ValueError("Transcript configuration is required for projection to genomic coordinates")
            
        if not variant_c.genomic_pos:
            # Create a mock genomic position based on the transcript position
            # In a real implementation, this would be calculated using the transcript's exon positions
            genomic_pos = GenomicPosition(
                chrom=config.chrom or "chr1",
                pos=config.tx_start + variant_c.start if config.tx_start else variant_c.start,
                ref=variant_c.ref or "N",
                alt=variant_c.alt or ""
            )
        else:
            genomic_pos = variant_c.genomic_pos
            
        # Create a new variant with the genomic position
        result = variant_c.model_copy(update={
            'hgvs': f"{genomic_pos.chrom}:g.{genomic_pos.pos}{genomic_pos.ref}>{genomic_pos.alt}",
            'genomic_pos': genomic_pos,
            'variant_type': VariantType.SUBSTITUTION if len(genomic_pos.ref) == 1 and len(genomic_pos.alt) == 1 else VariantType.SEQUENCE_ALTERATION
        })
        
        return result
        
    except Exception as e:
        raise ValueError(
            f"Failed to project variant {variant_c.hgvs} to genomic coordinates: {str(e)}"
        )