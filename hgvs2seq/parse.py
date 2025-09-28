"""
Handles parsing and normalization of HGVS variant strings.
"""
import hgvs.parser
import logging
from typing import List
from .models import VariantIn, VariantNorm, TranscriptConfig
from .data_provider import get_am, get_hdp

_logger = logging.getLogger(__name__)
_parser = hgvs.parser.Parser()

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
    # Get data providers
    am = get_am()
    hdp = get_hdp()

    normalized_variants = []
    transcript_id = config.transcript_id

    _logger.info(f"Parsing and normalizing {len(variants_in)} variants for {transcript_id}")

    for var_in in variants_in:
        try:
            # 1. Parse the HGVS string
            variant = _parser.parse_hgvs_variant(var_in.hgvs)
            _logger.debug(f"Parsed {var_in.hgvs} -> {variant}")

            # 2. Project to transcript coordinates (c.) if not already
            if variant.ac != transcript_id or variant.type != 'c':
                _logger.info(f"Projecting variant {variant} to {transcript_id}")
                variant_c = am.g_to_c(variant, transcript_id)
                _logger.info(f"Projected to {variant_c}")
            else:
                variant_c = variant

            # Ensure the variant is fully normalized by the data provider
            variant_norm = hdp.normalize_variant(variant_c)
            _logger.debug(f"Normalized {variant_c} -> {variant_norm}")

            # 3. Extract details into our VariantNorm model
            pos_edit = variant_norm.posedit
            edit = pos_edit.edit

            # Determine variant kind, distinguishing sub from mnv
            kind = edit.type
            if kind == 'sub':
                # A 'sub' from hgvs can be single-nucleotide or multi-nucleotide.
                # We classify multi-nucleotide subs as 'mnv'.
                if pos_edit.pos.start.base != pos_edit.pos.end.base:
                    kind = 'mnv'
                elif hasattr(edit, 'ref') and edit.ref and len(edit.ref) > 1:
                    kind = 'mnv'

            # Ensure the kind is one we explicitly support in our model
            from .models import VariantNorm
            valid_kinds = VariantNorm.model_fields['kind'].annotation.__args__
            if kind not in valid_kinds:
                _logger.error(f"Unsupported HGVS edit type '{kind}' for {var_in.hgvs}.")
                raise ValueError(f"Unsupported HGVS edit type: {kind}")

            # Start/end coordinates and offsets
            c_start = pos_edit.pos.start.base
            c_start_offset = pos_edit.pos.start.offset if hasattr(pos_edit.pos.start, 'offset') else 0
            c_end = pos_edit.pos.end.base
            c_end_offset = pos_edit.pos.end.offset if hasattr(pos_edit.pos.end, 'offset') else 0

            # Alt sequence
            alt = edit.alt if hasattr(edit, 'alt') and edit.alt is not None else ""

            norm_obj = VariantNorm(
                hgvs_c=str(variant_norm),
                kind=kind,
                c_start=c_start,
                c_start_offset=c_start_offset,
                c_end=c_end,
                c_end_offset=c_end_offset,
                alt=alt,
                meta={"original_hgvs": var_in.hgvs},
                phase_group=var_in.phase_group
            )
            normalized_variants.append(norm_obj)

        except Exception as e:
            _logger.error(f"Failed to parse or normalize variant '{var_in.hgvs}': {e}")
            raise ValueError(
                f"Failed to process variant '{var_in.hgvs}'. Ensure it is a valid "
                f"HGVS string and projectable to transcript {transcript_id}."
            ) from e

    _logger.info(f"Successfully normalized {len(normalized_variants)} variants.")
    return normalized_variants