import json
import logging
from pathlib import Path
from typing import TextIO

import click

from hgvs2seq import apply_variants, load_config, VariantIn
from hgvs2seq.config import TranscriptConfig

# Set up basic logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def parse_variants_file(file: TextIO) -> list[VariantIn]:
    """Parses a file with one HGVS string per line."""
    variants = []
    for line in file:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Rudimentary phase parsing from comments, e.g., #phase=1
        phase_group = None
        if "#" in line:
            comment = line.split("#", 1)[1]
            if "phase=" in comment:
                try:
                    phase_group = int(comment.split("phase=")[1].strip())
                except (ValueError, IndexError):
                    logger.warning(f"Could not parse phase from comment: '{comment}'")
            line = line.split("#", 1)[0].strip()

        variants.append(VariantIn(hgvs=line, phase_group=phase_group))
    return variants

@click.command()
@click.option(
    "--transcript",
    required=True,
    help="Transcript ID (e.g., NM_000551.3). A corresponding config JSON file is expected.",
)
@click.option(
    "--assembly",
    required=True,
    help="Reference assembly (e.g., GRCh38).",
)
@click.option(
    "--variants",
    "variants_file",
    type=click.File("r"),
    required=True,
    help="Path to a file containing HGVS strings, one per line.",
)
@click.option(
    "--out",
    "out_json_file",
    type=click.File("w"),
    help="Path to write the output JSON SequenceBundle.",
)
@click.option(
    "--fasta",
    "out_fasta_file",
    type=click.File("w"),
    help="Path to write the output FASTA sequences (cDNA and protein).",
)
@click.option(
    "--policy",
    type=click.Choice(["order_by_pos", "reject_overlaps"], case_sensitive=False),
    default="order_by_pos",
    show_default=True,
    help="Policy for handling overlapping variants.",
)
@click.option(
    "--emit",
    default="cdna,protein",
    show_default=True,
    help="Comma-separated list of sequences to output (e.g., 'cdna,protein').",
)
@click.option(
    "--annotate-splicing/--no-annotate-splicing",
    default=False,
    show_default=True,
    help="Enable SpliceAI annotation.",
)
@click.option(
    "--reference-path",
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the reference genome FASTA file (required for SpliceAI).",
)
@click.option(
    "--spliceai-distance",
    type=int,
    default=500,
    show_default=True,
    help="SpliceAI: maximum distance between variant and splice site.",
)
@click.option(
    "--annotate-nmd/--no-annotate-nmd",
    default=False,
    show_default=True,
    help="Enable NMD prediction.",
)
@click.option(
    "--nmd-threshold",
    type=int,
    default=50,
    show_default=True,
    help="NMD: distance from PTC to last exon-exon junction to trigger NMD.",
)
def main(
    transcript: str,
    assembly: str,
    variants_file: TextIO,
    out_json_file: TextIO | None,
    out_fasta_file: TextIO | None,
    policy: str,
    emit: str,
    annotate_splicing: bool,
    reference_path: str | None,
    spliceai_distance: int,
    annotate_nmd: bool,
    nmd_threshold: int,
):
    """
    Applies HGVS variants to a transcript and predicts the consequences.
    """
    logger.info(f"Processing transcript {transcript} on assembly {assembly}")

    # --- Configuration and Input Parsing ---
    try:
        # Assuming load_config can find a file like `NM_000551.3.json` or fetch it.
        config_path = f"{transcript}.json"
        if not Path(config_path).exists():
             logger.warning(f"Config file {config_path} not found. `load_config` will try to fetch from UTA.")
        cfg = load_config(config_path)
        # Override assembly if provided on CLI
        cfg.assembly = assembly
    except Exception as e:
        logger.error(f"Failed to load transcript configuration for {transcript}: {e}")
        raise click.Abort()

    variants = parse_variants_file(variants_file)
    logger.info(f"Loaded {len(variants)} variants from {variants_file.name}")

    outputs_set = set(emit.split(','))

    # --- Parameter Validation and Preparation ---
    if annotate_splicing and not reference_path:
        logger.error("--reference-path is required when --annotate-splicing is enabled.")
        raise click.Abort()

    spliceai_params = {
        "reference_path": reference_path,
        "distance": spliceai_distance,
    }

    nmd_params = {
        "nmd_threshold": nmd_threshold,
    }

    # --- Core Logic ---
    try:
        bundle = apply_variants(
            cfg=cfg,
            variants=variants,
            policy=policy,
            outputs=outputs_set,
            annotate_splicing=annotate_splicing,
            annotate_nmd=annotate_nmd,
            spliceai_params=spliceai_params,
            nmd_params=nmd_params,
        )
    except Exception as e:
        logger.error(f"An error occurred during variant application: {e}")
        raise click.Abort()

    # --- Output Generation ---
    if out_json_file:
        logger.info(f"Writing JSON output to {out_json_file.name}")
        out_json_file.write(bundle.model_dump_json(indent=2))

    if out_fasta_file:
        logger.info(f"Writing FASTA output to {out_fasta_file.name}")
        # Write reference sequences
        if "cdna" in outputs_set:
            out_fasta_file.write(f">cdna_ref|{cfg.transcript_id}\n")
            out_fasta_file.write(bundle.cdna_ref + "\n")
        if "protein" in outputs_set and bundle.protein_ref:
            out_fasta_file.write(f">protein_ref|{cfg.transcript_id}\n")
            out_fasta_file.write(bundle.protein_ref + "\n")

        # Write edited sequences per haplotype
        for i, (cdna_edited, prot_edited) in enumerate(zip(bundle.cdna_edited, bundle.protein_edited)):
            hap_id = f"haplotype_{i+1}"
            if "cdna" in outputs_set and cdna_edited:
                out_fasta_file.write(f">cdna_edited|{cfg.transcript_id}|{hap_id}\n")
                out_fasta_file.write(cdna_edited + "\n")
            if "protein" in outputs_set and prot_edited:
                out_fasta_file.write(f">protein_edited|{cfg.transcript_id}|{hap_id}\n")
                out_fasta_file.write(prot_edited + "\n")

    logger.info("Processing complete.")

if __name__ == "__main__":
    main()