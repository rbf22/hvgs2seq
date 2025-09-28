from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..parse import VariantNorm

def apply_single_variant(ref_seq: str, variant: "VariantNorm") -> str:
    """
    Applies a single normalized variant to a reference cDNA sequence.

    Args:
        ref_seq: The reference cDNA sequence.
        variant: A VariantNorm object representing the change.

    Returns:
        The modified cDNA sequence.

    Raises:
        ValueError: If the variant type is unknown or if the reference sequence
                    is inconsistent with the variant's position.
    """
    # HGVS coordinates are 1-based, Python strings are 0-based.
    # start_idx is inclusive, end_idx is exclusive for slicing.
    start_idx = variant.c_start - 1
    end_idx = variant.c_end

    # Basic bounds checking
    if not (0 <= start_idx <= end_idx <= len(ref_seq)):
        raise ValueError(
            f"Variant {variant.hgvs_c} coordinates ({variant.c_start}, {variant.c_end}) "
            f"are out of bounds for sequence of length {len(ref_seq)}."
        )

    # Dispatch based on variant kind
    kind = variant.kind
    alt = variant.alt

    if kind == "sub":
        if alt is None:
            raise ValueError(f"Substitution variant {variant.hgvs_c} is missing 'alt' sequence.")
        # Substitution of a single base
        return ref_seq[:start_idx] + alt + ref_seq[end_idx:]

    elif kind == "del":
        # Deletion of the sequence between start and end
        return ref_seq[:start_idx] + ref_seq[end_idx:]

    elif kind == "ins":
        if alt is None:
            raise ValueError(f"Insertion variant {variant.hgvs_c} is missing 'alt' sequence.")
        # Insertion happens between two bases. HGVS c.1_2insA means insert after base 1.
        # c_start is the base 5' to the insertion site, so we insert at its 0-based index.
        # In this case, c_start=1, so index is 1.
        # Let's re-verify: hgvs-python for c.1_2insA gives start=1, end=2.
        # The insertion point is between them. In 0-based, that's after index 0, at index 1.
        # So we insert at `start_idx + 1`.
        # No, c.1_2insA means between base 1 and 2. So after index 0, and before index 1.
        # The correct insertion point is `start_idx + 1`, but `hgvs-python` might give c_start as the insertion point.
        # Let's assume c_start is the position *before* the insertion.
        # `c.123_124insT` means insertion is at index 123. c_start=123.
        # So, the insertion point is `start_idx`.
        # Let's stick to the simplest interpretation: for c.1_2insA, c_start=1, c_end=2.
        # Insertion point is after c_start. So after index 0. At index 1.
        # `ref_seq[:1] + alt + ref_seq[1:]`. This seems right.
        # The insertion point is `c_end - 1`
        insertion_idx = variant.c_end - 1
        return ref_seq[:insertion_idx] + alt + ref_seq[insertion_idx:]


    elif kind in ("delins", "mnv"):
        if alt is None:
            raise ValueError(f"Variant {variant.hgvs_c} is missing 'alt' sequence.")
        # Deletion of the specified range, followed by an insertion at the start of that range
        return ref_seq[:start_idx] + alt + ref_seq[end_idx:]

    elif kind == "dup":
        # Duplication of the sequence in the given range
        duplicated_seq = ref_seq[start_idx:end_idx]
        # The duplication is inserted immediately after the original sequence
        return ref_seq[:end_idx] + duplicated_seq + ref_seq[end_idx:]

    else:
        # Currently, 'inv' (inversion) is not supported.
        raise ValueError(f"Unsupported variant kind: '{kind}'")