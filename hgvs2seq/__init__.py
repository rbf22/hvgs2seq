from typing import List, Literal, Set, Dict, Any, Optional

from .config import TranscriptConfig, load_config
from .models import SequenceBundle, EditPlan
from .parse import VariantIn, VariantNorm, normalize_variant
from .refseq import get_reference_sequence
from .apply.batch import apply_edit_plan
from .consequence.cds import get_protein_sequence_and_consequences
from .splicing.spliceai import annotate_spliceai
from .consequence.nmd import check_nmd

__all__ = ["apply_variants", "load_config", "TranscriptConfig", "VariantIn", "SequenceBundle"]

def apply_variants(
    cfg: TranscriptConfig,
    variants: List[VariantIn],
    *,
    policy: Literal["order_by_pos", "reject_overlaps"] = "order_by_pos",
    outputs: Set[str] = {"cdna", "protein"},
    annotate_splicing: bool = False,
    annotate_nmd: bool = False,
    spliceai_params: Optional[Dict[str, Any]] = None,
) -> SequenceBundle:
    """
    The main entry point for applying a list of HGVS variants to a transcript.

    This function orchestrates the entire process:
    1. Fetches the reference sequence.
    2. Parses and normalizes input variants.
    3. Groups variants into haplotypes.
    4. Applies variants to the reference sequence.
    5. Computes the consequences on the protein sequence.
    6. Optionally, annotates for splicing effects and NMD likelihood.

    Args:
        cfg: The transcript configuration.
        variants: A list of input variants to apply.
        policy: The policy for handling overlapping variants within a haplotype.
        outputs: A set specifying which sequences to include in the output (e.g., "cdna", "protein").
        annotate_splicing: If True, run SpliceAI annotations.
        annotate_nmd: If True, run NMD checks.
        spliceai_params: Optional parameters for the SpliceAI annotation function.
                         Must include 'reference_path' if annotate_splicing is True.

    Returns:
        A SequenceBundle containing the reference and edited sequences, annotations, and provenance.
    """
    # 1. Fetch reference sequence
    ref_cdna = get_reference_sequence(cfg.transcript_id)

    # 2. Normalize all variants
    normalized_variants: List[VariantNorm] = []
    parsing_warnings = []
    for var_in in variants:
        try:
            norm_var = normalize_variant(var_in, cfg)
            norm_var.meta['phase_group'] = var_in.phase_group
            normalized_variants.append(norm_var)
        except ValueError as e:
            parsing_warnings.append(str(e))

    # 3. Group variants by phase group into haplotypes
    haplotype_map: Dict[int, List[VariantNorm]] = {}
    unphased = []
    for v in normalized_variants:
        phase_group = v.meta.get('phase_group')
        if phase_group is not None:
            if phase_group not in haplotype_map:
                haplotype_map[phase_group] = []
            haplotype_map[phase_group].append(v)
        else:
            unphased.append(v)

    final_haplotypes: List[List[VariantNorm]] = []
    if haplotype_map:
        for phase_group in sorted(haplotype_map.keys()):
            final_haplotypes.append(unphased + haplotype_map[phase_group])
    else:
        final_haplotypes.append(unphased)

    # 4. Create an EditPlan
    edit_plan = EditPlan(
        haplotypes=final_haplotypes,
        policy=policy,
        warnings=parsing_warnings,
    )

    # 5. Apply variants to get edited cDNA sequences
    edited_cdn_sequences = apply_edit_plan(ref_cdna, edit_plan)

    # 6. Compute protein consequences and other annotations
    protein_ref: Optional[str] = None
    protein_edited_list: List[Optional[str]] = []
    annotations: Dict[str, Any] = {"haplotypes": []}

    if "protein" in outputs:
        protein_ref, _, _ = get_protein_sequence_and_consequences(ref_cdna, ref_cdna, cfg)

    splicing_annotations = {}
    if annotate_splicing:
        splicing_annotations = annotate_spliceai(
            variants=normalized_variants,
            config=cfg,
            **(spliceai_params or {})
        )

    for i, edited_cdna in enumerate(edited_cdn_sequences):
        hap_ann: Dict[str, Any] = {"haplotype_index": i + 1}
        prot_edited: Optional[str] = None

        if edited_cdna:
            if "protein" in outputs:
                prot_ref_calc, prot_edited, cons = get_protein_sequence_and_consequences(
                    ref_cdna, edited_cdna, cfg
                )
                if protein_ref is None: protein_ref = prot_ref_calc
                hap_ann["protein_consequence"] = cons

            if annotate_splicing:
                haplotype_variants = edit_plan.haplotypes[i]
                hap_splicing_ann = {
                    v.hgvs_c: splicing_annotations.get(v.hgvs_c, {"status": "not_found"})
                    for v in haplotype_variants
                }
                hap_ann["splicing"] = hap_splicing_ann

            if annotate_nmd:
                hap_ann["nmd"] = check_nmd(edited_cdna, protein_ref, cfg)

        protein_edited_list.append(prot_edited)
        annotations["haplotypes"].append(hap_ann)

    # 7. Assemble the final SequenceBundle
    bundle = SequenceBundle(
        cdna_ref=ref_cdna,
        cdna_edited=edited_cdn_sequences,
        mrna_ref=ref_cdna,
        mrna_edited=edited_cdn_sequences,
        protein_ref=protein_ref,
        protein_edited=protein_edited_list,
        annotations=annotations,
        provenance={
            "transcript_id": cfg.transcript_id,
            "assembly": cfg.assembly,
            "variants_input": [v.model_dump() for v in variants],
            "variants_normalized": [v.model_dump() for v in normalized_variants],
            "edit_plan_warnings": edit_plan.warnings,
        },
    )

    return bundle