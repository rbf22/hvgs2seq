import subprocess
import tempfile
import logging
from typing import List, Dict, Any, cast

import hgvs.parser
import hgvs.sequencevariant
import hgvs.posedit
import hgvs.location

from ..models import VariantNorm
from ..config import TranscriptConfig
from ..project import project_c_to_g

logger = logging.getLogger(__name__)
_hp = hgvs.parser.Parser()


def _parse_spliceai_info(info_str: str) -> Dict[str, Any]:
    """
    Parses the SpliceAI annotation from the VCF INFO field.
    The format is: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    """
    try:
        # The annotation is typically prefixed with 'SpliceAI='
        spliceai_field = next(s for s in info_str.split(';') if s.startswith('SpliceAI='))
        _, data = spliceai_field.split('=', 1)

        # There can be multiple predictions for a variant (e.g., if it hits multiple genes).
        # We'll take the first one for simplicity here. A more robust implementation might handle all.
        first_prediction = data.split(',')[0]
        parts = first_prediction.split('|')

        return {
            "allele": parts[0],
            "symbol": parts[1],
            "ds_ag": float(parts[2]), # delta score acceptor gain
            "ds_al": float(parts[3]), # delta score acceptor loss
            "ds_dg": float(parts[4]), # delta score donor gain
            "ds_dl": float(parts[5]), # delta score donor loss
            "dp_ag": int(parts[6]),   # delta position acceptor gain
            "dp_al": int(parts[7]),   # delta position acceptor loss
            "dp_dg": int(parts[8]),   # delta position donor gain
            "dp_dl": int(parts[9]),   # delta position donor loss
        }
    except (StopIteration, IndexError, ValueError) as e:
        logger.warning(f"Could not parse SpliceAI info field: '{info_str}'. Error: {e}")
        return {"status": "parsing_failed", "raw_info": info_str}


def annotate_spliceai(
    variants: List[VariantNorm],
    config: "TranscriptConfig",
    *,
    reference_path: str,
    distance: int = 500,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    Annotates variants with SpliceAI predictions by calling the command-line tool.

    Args:
        variants: A list of normalized variants to annotate.
        config: The transcript configuration.
        reference_path: The file path to the reference genome FASTA file (e.g., hg38.fa).
                        This is a required argument.
        distance: Maximum distance between the variant and splice site for annotation.
        **kwargs: Other keyword arguments (currently unused).

    Returns:
        A dictionary of annotations for each variant.
    """
    if not reference_path:
        raise ValueError("A valid 'reference_path' to a FASTA file is required for SpliceAI.")

    annotations = {}
    if config.assembly.lower() not in ["grch37", "grch38", "hg19", "hg38"]:
         raise ValueError(f"Unsupported assembly for SpliceAI: {config.assembly}. Only grch37/hg19 and grch38/hg38 are supported by the bundled annotations.")

    assembly_arg = "grch38" if config.assembly.lower() in ["grch38", "hg38"] else "grch37"

    # Create VCF content in-memory
    vcf_header = "##fileformat=VCFv4.2\n"
    vcf_header += f"##reference={config.assembly}\n"
    vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

    vcf_records = []
    variants_to_process = []
    for variant in variants:
        try:
            # Project c. variant to g. coordinate
            var_c = _hp.parse_hgvs_variant(variant.hgvs_c)
            var_g = project_c_to_g(var_c)

            pos_g = cast(hgvs.location.SimplePosition, var_g.posedit.pos.start)
            edit_g = var_g.posedit.edit

            chrom = var_g.ac
            pos = pos_g.base
            ref = edit_g.ref
            alt = edit_g.alt

            # Use the c. HGVS string as the ID to map results back
            vcf_records.append(f"{chrom}\t{pos}\t{variant.hgvs_c}\t{ref}\t{alt}\t.\tPASS\t.\n")
            variants_to_process.append(variant)
        except (ValueError, hgvs.exceptions.HGVSError) as e:
            logger.error(f"Could not project variant {variant.hgvs_c} for SpliceAI: {e}")
            annotations[variant.hgvs_c] = {"status": "projection_failed", "error": str(e)}

    if not vcf_records:
        logger.warning("No variants could be processed for SpliceAI.")
        return annotations

    vcf_content = vcf_header + "".join(vcf_records)

    with tempfile.NamedTemporaryFile(mode='w', delete=True, suffix=".vcf") as input_vcf, \
         tempfile.NamedTemporaryFile(mode='r', delete=True, suffix=".vcf") as output_vcf:

        input_vcf.write(vcf_content)
        input_vcf.flush()

        command = [
            "spliceai",
            "-I", input_vcf.name,
            "-O", output_vcf.name,
            "-R", reference_path, # This needs to be a path to the reference fasta file
            "-A", assembly_arg,
            "-D", str(distance)
        ]

        try:
            # SpliceAI can be noisy, so we capture stderr
            logger.info(f"Running SpliceAI with command: {' '.join(command)}")
            result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=300)
            logger.debug(f"SpliceAI stderr:\n{result.stderr}")

            # Parse the output VCF
            output_vcf.seek(0)
            for line in output_vcf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                var_id = fields[2]
                info = fields[7]
                annotations[var_id] = _parse_spliceai_info(info)

        except FileNotFoundError:
            msg = "The 'spliceai' command was not found. Is it installed and in the system's PATH?"
            logger.error(msg)
            raise RuntimeError(msg)
        except subprocess.CalledProcessError as e:
            msg = f"SpliceAI execution failed with exit code {e.returncode}.\n"
            msg += f"Stderr: {e.stderr}\n"
            msg += "This may be due to a missing reference genome file or incorrect VCF input."
            logger.error(msg)
            # Annotate all processed variants with the failure
            for var in variants_to_process:
                annotations[var.hgvs_c] = {"status": "execution_failed", "error": msg}
        except Exception as e:
            logger.error(f"An unexpected error occurred during SpliceAI annotation: {e}")
            for var in variants_to_process:
                annotations[var.hgvs_c] = {"status": "unexpected_error", "error": str(e)}

    return annotations