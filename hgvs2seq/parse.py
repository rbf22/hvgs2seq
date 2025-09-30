"""
Handles parsing and normalization of HGVS variant strings.
"""
import logging
from typing import List, Optional
from .models import VariantIn, VariantNorm, TranscriptConfig
from .data_provider import get_transcript_data, get_genome_sequence
from .hgvs_compat import (
    parser, 
    SequenceVariant, 
    HGVSParseError, 
    HGVSDataNotAvailableError,
    Interval,  # Add Interval to imports
    PosEdit,   # Add PosEdit to imports
    SimplePosition  # Add SimplePosition to imports
)

_logger = logging.getLogger(__name__)

# Use the parser from our compatibility layer
_parser = parser

def parse_and_normalize_variants(
    variants_in: List[VariantIn],
    config: TranscriptConfig
) -> List[VariantNorm]:
    """
    Parses a list of HGVS strings, projects them to the specified transcript,
    and normalizes them into a structured format.

    Args:
        variants_in: A list of VariantIn objects containing HGVS strings.
        config: The TranscriptConfig for the target transcript.

    Returns:
        A list of VariantNorm objects.

    Raises:
        ValueError: If a variant cannot be parsed or projected to the transcript.
    """
    normalized_variants = []
    transcript_id = config.transcript_id

    _logger.info(f"Parsing and normalizing {len(variants_in)} variants for {transcript_id}")

    for var_in in variants_in:
        try:
            # 1. Parse the HGVS string
            variant = _parser.parse_hgvs_variant(var_in.hgvs)
            _logger.debug(f"Parsed {var_in.hgvs} -> {variant}")

            # 2. Get transcript data
            transcript_data = get_transcript_data(transcript_id)
            if not transcript_data:
                raise HGVSDataNotAvailableError(
                    f"Transcript data not available for {transcript_id}"
                )

            # 3. Project to transcript coordinates (c.) if not already
            if variant.ac != transcript_id or variant.type != 'c':
                _logger.info(f"Projecting variant {variant} to {transcript_id}")
                
                # For pyhgvs, we'll need to implement the projection logic
                # This is a simplified version - in a real implementation, 
                # we would use pyhgvs functions to do the projection
                if variant.type == 'g':
                    # Project genomic to transcript coordinates
                    variant_c = _project_genomic_to_transcript(variant, transcript_data)
                elif variant.type == 'n':
                    # Project non-coding to coding (if applicable)
                    variant_c = _project_noncoding_to_coding(variant, transcript_data)
                elif variant.type == 'c':
                    # Convert between transcripts
                    variant_c = _convert_between_transcripts(variant, transcript_data)
                else:
                    raise ValueError(f"Unsupported variant type: {variant.type}")
                    
                _logger.info(f"Projected to {variant_c}")
            else:
                variant_c = variant

            # Normalize the variant (simplified for now)
            variant_norm = _normalize_variant(variant_c, transcript_data)
            _logger.debug(f"Normalized {variant_c} -> {variant_norm}")

            # Create the normalized variant
            # Get the variant type from the edit object
            variant_type = variant_norm.posedit.edit.type if hasattr(variant_norm.posedit, 'edit') and hasattr(variant_norm.posedit.edit, 'type') else 'alt'
            
            normalized = VariantNorm(
                hgvs=var_in.hgvs,
                transcript_id=variant_norm.ac,
                variant_type=variant_type,  # This will be converted to VariantType by the model
                start=variant_norm.posedit.pos.start.base,
                end=variant_norm.posedit.pos.end.base,
                ref=variant_norm.posedit.edit.ref if hasattr(variant_norm.posedit.edit, 'ref') else "",
                alt=variant_norm.posedit.edit.alt if hasattr(variant_norm.posedit.edit, 'alt') else "",
                hgvs_c=str(variant_norm),
                hgvs_p="",  # Will be filled in later if needed
                consequence="",  # Will be determined later
                is_canonical=True,  # Default to true, can be updated later
                phase_group=var_in.phase_group  # Pass through the phase group
            )
            normalized_variants.append(normalized)

        except HGVSParseError as e:
            _logger.error(f"Failed to parse HGVS string '{var_in.hgvs}': {e}")
            raise ValueError(f"Invalid HGVS string: {var_in.hgvs}") from e
        except HGVSDataNotAvailableError as e:
            _logger.error(f"Data not available for variant '{var_in.hgvs}': {e}")
            raise ValueError(
                f"Could not process variant {var_in.hgvs}: Required data not available"
            ) from e
        except Exception as e:
            _logger.error(f"Unexpected error processing variant '{var_in.hgvs}': {e}")
            raise ValueError(
                f"Failed to process variant {var_in.hgvs}: {str(e)}"
            ) from e

    return normalized_variants


def _project_genomic_to_transcript(
    variant: SequenceVariant, 
    transcript_data: dict
) -> SequenceVariant:
    """
    Project a genomic variant to transcript coordinates.
    
    Args:
        variant: The genomic variant to project
        transcript_data: Transcript data from the data provider
        
    Returns:
        A new SequenceVariant in transcript coordinates
    """
    # This is a simplified implementation
    # In a real implementation, we would use pyhgvs functions to do the projection
    # based on the transcript's exon structure and sequence
    
    # For now, we'll just create a new variant with the transcript ID
    # and the same position (this is not correct, just a placeholder)
    return SequenceVariant(
        ac=transcript_data['transcript_id'],
        type='c',
        posedit=variant.posedit
    )


def _project_noncoding_to_coding(
    variant: SequenceVariant, 
    transcript_data: dict
) -> SequenceVariant:
    """
    Project a non-coding variant to coding coordinates if possible.
    
    Args:
        variant: The non-coding variant to project
        transcript_data: Transcript data from the data provider
        
    Returns:
        A new SequenceVariant in coding coordinates
    """
    # This is a simplified implementation
    # In a real implementation, we would use the transcript's CDS information
    # to convert from non-coding to coding coordinates
    
    return SequenceVariant(
        ac=variant.ac,
        type='c',
        posedit=variant.posedit
    )


def _convert_between_transcripts(
    variant: SequenceVariant, 
    target_transcript_data: dict
) -> SequenceVariant:
    """
    Convert a variant from one transcript to another.
    
    Args:
        variant: The variant to convert
        target_transcript_data: Target transcript data
        
    Returns:
        A new SequenceVariant in the target transcript's coordinates
    """
    # This is a simplified implementation
    # In a real implementation, we would use the transcript sequences
    # to map between different transcripts
    
    return SequenceVariant(
        ac=target_transcript_data['transcript_id'],
        type='c',
        posedit=variant.posedit
    )


def _normalize_variant(
    variant: SequenceVariant,
    transcript_data: dict
) -> SequenceVariant:
    """
    Normalize a variant to its most concise representation.
    
    Args:
        variant: The variant to normalize
        transcript_data: Transcript data for context
        
    Returns:
        A new normalized SequenceVariant
    """
    # This is a simplified implementation
    # In a real implementation, we would normalize the variant
    # by left-aligning and trimming common bases
    
    # For now, just return a copy of the variant
    # with some basic normalization
    pos = variant.posedit.pos
    
    # Create a normalized position
    normalized_pos = Interval(
        start=SimplePosition(base=min(pos.start.base, pos.end.base)),
        end=SimplePosition(base=max(pos.start.base, pos.end.base)),
        uncertain=pos.uncertain
    )
    
    # Create a new variant with the normalized position
    return SequenceVariant(
        ac=variant.ac,
        type=variant.type,
        posedit=PosEdit(pos=normalized_pos, edit=variant.posedit.edit)
    )