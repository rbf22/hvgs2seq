from typing import TYPE_CHECKING
import hgvs
import hgvs.dataproviders.uta
import hgvs.variantmapper
from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError

# Use a forward reference for TranscriptConfig to avoid circular imports at runtime
if TYPE_CHECKING:
    from .config import TranscriptConfig

# Connect to the public UTA instance.
# In a production application, this might be made more configurable,
# or could use a local SeqRepo instance for better performance and reproducibility.
try:
    hdp = hgvs.dataproviders.uta.connect()
    # The VariantMapper is the core tool for projecting variants between coordinates.
    mapper = hgvs.variantmapper.VariantMapper(hdp)
except Exception as e:
    # If the connection fails on startup, we can't proceed.
    # We'll raise a more informative error here.
    raise RuntimeError(
        "Failed to connect to HGVS data provider (UTA). "
        "Please check network connection and UTA service status. "
        f"Original error: {e}"
    )


def project_variant(
    variant: hgvs.sequencevariant.SequenceVariant, config: "TranscriptConfig"
) -> hgvs.sequencevariant.SequenceVariant:
    """
    Projects a variant to the specified transcript's cDNA (c.) coordinates.

    This function takes any valid HGVS variant and attempts to convert it to the
    cDNA coordinates of the transcript defined in the configuration.

    :param variant: An hgvs SequenceVariant object.
    :param config: The TranscriptConfig for the target transcript.
    :return: A new SequenceVariant object in c. coordinates.
    :raises ValueError: If the variant cannot be projected to the target transcript.
    """
    target_ac = config.transcript_id

    # If the variant is already in the correct form, we're done.
    if variant.ac == target_ac and variant.type == "c":
        return variant

    try:
        # Dispatch to the correct mapper method based on the variant type.
        if variant.type == "g":
            var_c = mapper.g_to_c(variant, target_ac)
        elif variant.type == "n":
            var_c = mapper.n_to_c(variant, target_ac)
        elif variant.type == "c":
            # This handles projections between different transcripts.
            var_c = mapper.c_to_c(variant, target_ac)
        else:
            raise ValueError(f"Unsupported variant type for projection: '{variant.type}'")
        return var_c
    except HGVSDataNotAvailableError as e:
        # This error occurs if UTA doesn't have the necessary data to perform the mapping.
        raise ValueError(
            f"Could not project variant {variant} to {target_ac}. "
            f"Data not available from the data provider: {e}"
        )
    except HGVSError as e:
        # This is a general catch-all for other mapping errors from the hgvs library.
        raise ValueError(
            f"Failed to project variant {variant} to {target_ac} due to a mapping error: {e}"
        )


def project_c_to_g(
    variant_c: hgvs.sequencevariant.SequenceVariant,
) -> hgvs.sequencevariant.SequenceVariant:
    """
    Projects a cDNA (c.) variant to genomic (g.) coordinates.

    This function requires that the input variant is in c. coordinates and
    is associated with a transcript that can be found by the data provider.

    :param variant_c: An hgvs SequenceVariant object in c. coordinates.
    :return: A new SequenceVariant object in g. coordinates.
    :raises ValueError: If the variant cannot be projected.
    """
    if variant_c.type != "c":
        raise ValueError(f"Input variant must be of type 'c' for c_to_g projection, but got '{variant_c.type}'.")

    try:
        # The mapper uses the transcript information from the variant's accession (e.g., NM_...)
        # to determine the correct genomic context.
        var_g = mapper.c_to_g(variant_c)
        return var_g
    except HGVSDataNotAvailableError as e:
        raise ValueError(
            f"Could not project variant {variant_c} to genomic coordinates. "
            f"Data not available from the data provider: {e}"
        )
    except HGVSError as e:
        raise ValueError(
            f"Failed to project variant {variant_c} to genomic coordinates due to a mapping error: {e}"
        )