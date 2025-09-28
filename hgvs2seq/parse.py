from pydantic import BaseModel
from typing import Literal, Optional, TYPE_CHECKING, cast

import hgvs.parser
import hgvs.edit
import hgvs.posedit
import hgvs.location
from hgvs.exceptions import HGVSParseError

# Use a forward reference for TranscriptConfig to avoid circular imports at runtime
if TYPE_CHECKING:
    from .config import TranscriptConfig

from .project import project_variant


class VariantIn(BaseModel):
    """
    Represents an input variant with an HGVS string and optional phasing information.
    """
    hgvs: str
    phase_group: Optional[int] = None


class VariantNorm(BaseModel):
    """
    Represents a normalized variant in cDNA (c.) coordinates.
    """
    hgvs_c: str
    kind: Literal["sub", "del", "ins", "delins", "dup", "inv", "mnv"]
    c_start: int
    c_end: int
    alt: Optional[str] = None  # Inserted/alternate bases, normalized. None for deletions.
    meta: dict  # Original string, notes, etc.


# Initialize the hgvs parser. It's stateless and can be shared.
_hp = hgvs.parser.Parser()


def _get_kind_from_edit(edit: hgvs.edit.Edit) -> Literal["sub", "del", "ins", "delins", "dup", "inv", "mnv"]:
    """Maps an hgvs Edit object to our internal 'kind' string."""
    if isinstance(edit, hgvs.edit.Dup):
        return "dup"
    elif isinstance(edit, hgvs.edit.Inv):
        return "inv"
    elif isinstance(edit, hgvs.edit.NARefAlt):
        # The `type` property on NARefAlt gives us one of: "sub", "del", "ins", "delins", "identity"
        edit_type = edit.type

        if edit_type == "identity":
            raise ValueError("Identity variants (e.g., 'c.1A=A') are not supported.")

        # Differentiate between a generic delins and a multi-nucleotide variant (MNV),
        # which we define as a delins of the same length (e.g., c.1_2delAGinsTC)
        if edit_type == "delins":
            if edit.ref and edit.alt and len(edit.ref) == len(edit.alt) and len(edit.ref) > 1:
                return "mnv"

        # The type property value ("sub", "del", "ins", "delins") directly maps to our kind.
        return edit_type

    raise TypeError(f"Unknown or unsupported edit type: {type(edit)}")


def normalize_variant(variant_in: VariantIn, config: "TranscriptConfig") -> VariantNorm:
    """
    Parses an HGVS string, projects it to the target transcript's c. coordinates,
    and returns a structured VariantNorm object.
    """
    try:
        raw_variant = _hp.parse_hgvs_variant(variant_in.hgvs)
    except HGVSParseError as e:
        raise ValueError(f"Failed to parse HGVS string '{variant_in.hgvs}': {e}")

    # Project the variant to the target transcript's c. coordinates
    var_c = project_variant(raw_variant, config)

    # Extract details from the normalized hgvs object
    posedit = cast(hgvs.posedit.PosEdit, var_c.posedit)
    edit = posedit.edit
    pos = posedit.pos

    kind = _get_kind_from_edit(edit)
    alt = getattr(edit, 'alt', None)

    # HGVS locations are 1-based
    start_pos = cast(hgvs.location.SimplePosition, pos.start)
    end_pos = cast(hgvs.location.SimplePosition, pos.end)

    return VariantNorm(
        hgvs_c=str(var_c),
        kind=kind,
        c_start=start_pos.base,
        c_end=end_pos.base,
        alt=alt,
        meta={"original_hgvs": variant_in.hgvs},
    )