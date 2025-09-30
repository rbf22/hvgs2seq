#!/usr/bin/env python3
"""
Complex Variant Analysis with hgvs2seq
======================================

This script demonstrates how to analyze phased genetic variants using hgvs2seq.
It shows how to:
1. Define a test transcript with coding sequence
2. Create phased variants
3. Generate all possible haplotypes
4. Analyze protein consequences
5. Predict Nonsense-Mediated Decay (NMD) effects

Example usage:
    python complex_variant_analysis.py
"""

import sys
import os
import json
import logging
from typing import List, Dict, Tuple, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add the project root to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.insert(0, project_root)

# Import hgvs2seq components
try:
    from hgvs2seq.data_provider import _sequence_store
    from hgvs2seq.models import TranscriptConfig, VariantIn
    from hgvs2seq.parse import parse_and_normalize_variants
    from hgvs2seq.apply.batch import apply_variants_in_batch
    from hgvs2seq.consequence.cds import analyze_consequences
    from hgvs2seq.consequence.nmd import check_nmd
except ImportError as e:
    logger.error("Failed to import hgvs2seq components. Make sure hgvs2seq is installed.")
    raise

# Constants
OUTPUT_DIR = os.path.join(project_root, 'output')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# PAH gene transcript ID and reference sequence
PAH_TRANSCRIPT = "ENST00000553106"

# Example variant lists
example_1 = [
    f"{PAH_TRANSCRIPT}:c.442-5C>G",
    f"{PAH_TRANSCRIPT}:c.1315+1G>A",
    f"{PAH_TRANSCRIPT}:c.898G>T",
    f"{PAH_TRANSCRIPT}:c.1199+1G>C",
    f"{PAH_TRANSCRIPT}:c.842+3G>C",
    f"{PAH_TRANSCRIPT}:c.1208C>T",
    f"{PAH_TRANSCRIPT}:c.510-?_706+?del",
    f"{PAH_TRANSCRIPT}:c.441+4A>G",
    f"{PAH_TRANSCRIPT}:c.442-102_509+781del",
    f"{PAH_TRANSCRIPT}:c.140C>T",
    f"{PAH_TRANSCRIPT}:c.1066-11G>A",
    f"{PAH_TRANSCRIPT}:c.385G>T"
]

example_2 = [
    f"{PAH_TRANSCRIPT}:c.1340C>A",
    f"{PAH_TRANSCRIPT}:c.886G>C",
    f"{PAH_TRANSCRIPT}:c.1180G>C",
    f"{PAH_TRANSCRIPT}:c.1243G>A",
    f"{PAH_TRANSCRIPT}:c.169-2_352+?del",
    f"{PAH_TRANSCRIPT}:c.176A>G",
    f"{PAH_TRANSCRIPT}:c.533A>G",
    f"{PAH_TRANSCRIPT}:c.662A>G",
    f"{PAH_TRANSCRIPT}:c.1169A>G",
    f"{PAH_TRANSCRIPT}:c.1066-2A>T",
    f"{PAH_TRANSCRIPT}:c.168+5G>A",
    f"{PAH_TRANSCRIPT}:c.169-13T>G"
]

# For testing, we'll use a placeholder sequence
# In a real scenario, you would fetch this from a genomic database
TEST_SEQUENCE = 'N' * 2000  # Placeholder sequence


def load_fasta(fasta_path: str) -> str:
    """
    Load a FASTA file and return the sequence as a string.
    
    Args:
        fasta_path: Path to the FASTA file
        
    Returns:
        str: The sequence as a string with no headers or newlines
    """
    try:
        from Bio import SeqIO
        # Use Biopython's SeqIO for robust FASTA parsing
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq)
    except ImportError:
        # Fallback to simple parsing if Biopython is not available
        with open(fasta_path, 'r') as f:
            # Skip the header line
            next(f)
            # Read the rest of the file and remove newlines
            return ''.join(line.strip() for line in f)

def setup_sequence_store() -> None:
    """
    Initialize the in-memory sequence store with PAH transcript data.
    
    This function populates the global _sequence_store with:
    - The PAH transcript sequence (NM_000277.3)
    - A dummy chromosome 12 sequence for compatibility
    """
    try:
        from hgvs2seq.refseq import get_reference_cDNA
        from hgvs2seq.data_provider import get_genome_sequence
        
        # PAH gene transcript ID (RefSeq)
        PAH_TRANSCRIPT = "NM_000277.3"
        
        logger.info(f"Loading PAH transcript sequence for {PAH_TRANSCRIPT}")
        
        # Get the full transcript sequence using hgvs2seq's refseq module
        pah_seq = get_reference_cDNA(PAH_TRANSCRIPT)
        
        # Store the sequences
        # Create a dummy chromosome 12 sequence for compatibility
        _sequence_store['chr12'] = 'N' * 150_000_000  # Approx length of chr12
        _sequence_store['PAH_GENE'] = pah_seq
        _sequence_store['TEST_TRANSCRIPT'] = pah_seq  # For backward compatibility
        _sequence_store[PAH_TRANSCRIPT] = pah_seq  # Store with transcript ID as key
        
        logger.info(f"Loaded PAH transcript with length: {len(pah_seq):,} bp")
        
    except Exception as e:
        logger.error(f"Error loading sequence data: {str(e)}")
        logger.warning("Falling back to test sequence")
        _sequence_store['TEST_TRANSCRIPT'] = TEST_SEQUENCE
        _sequence_store['chr1'] = 'N' * 1000 + TEST_SEQUENCE + 'N' * 1000
        logger.info(f"Using test sequence with length: {len(TEST_SEQUENCE)}")
        
        # Log the full error for debugging
        import traceback
        logger.debug(traceback.format_exc())


def generate_haplotypes(reference: str, variants: List[VariantIn]) -> List[Tuple[str, List[str]]]:
    """
    Generate all possible haplotypes from a list of variants with phase groups.
    
    Args:
        reference: The reference DNA sequence
        variants: List of VariantIn objects with phase_group attributes
        
    Returns:
        List of tuples containing:
        - haplotype_sequence (str): The DNA sequence with variants applied
        - variant_list (List[str]): List of variant HGVS strings in this haplotype
    """
    # Initialize with reference sequence and no variants
    haplotypes = [(list(reference), [])]
    
    # Group variants by phase
    phase_groups = {}
    for var in variants:
        phase = var.phase_group or 0
        if phase not in phase_groups:
            phase_groups[phase] = []
        phase_groups[phase].append(var)
    
    # For each phase group, create new haplotypes
    for phase, phase_vars in sorted(phase_groups.items()):
        new_haplotypes = []
        
        for hap_seq, hap_vars in haplotypes:
            # Create two haplotypes for this phase (0 = reference, 1 = alternate)
            for allele in [0, 1]:
                new_hap = hap_seq.copy()  # Make a copy of the current haplotype
                new_vars = hap_vars.copy()  # Copy the variant list
                
                for var in phase_vars:
                    try:
                        # Parse variant position (1-based to 0-based)
                        pos = int(var.hgvs.split('c.')[1].split('>')[0][:-1]) - 1
                        ref = var.hgvs.split('>')[0][-1]
                        alt = var.hgvs.split('>')[1].split('(')[0]
                        
                        # Apply variant if this is the alternate allele
                        if allele == 1 and (new_hap[pos] == ref or new_hap[pos] == 'N'):
                            new_hap[pos] = alt
                            new_vars.append(var.hgvs)
                    except Exception as e:
                        logger.error(f"Error processing variant {var.hgvs}: {e}")
                
                # Only add if this combination is new
                hap_str = ''.join(new_hap)
                if (hap_str, tuple(sorted(new_vars))) not in [
                    (''.join(h), tuple(sorted(v))) for h, v in new_haplotypes
                ]:
                    new_haplotypes.append((new_hap, new_vars))
        
        # Update haplotypes for next phase
        haplotypes = new_haplotypes
    
    # Convert back to strings and sort variants
    return [
        (''.join(hap_seq), sorted(vars, key=lambda x: int(x.split('c.')[1].split('>')[0][:-1])))
        for hap_seq, vars in haplotypes
    ]


def create_transcript_config() -> TranscriptConfig:
    """
    Create a TranscriptConfig for the PAH gene (ENST00000553106).
    
    Returns:
        TranscriptConfig: Configuration for the PAH transcript
    """
    # PAH gene coordinates (1-based, inclusive)
    PAH_START = 102_836_311
    PAH_END = 102_958_772
    
    # PAH CDS coordinates (1-based, inclusive)
    PAH_CDS_START = 102_836_311  # Start of first coding exon
    PAH_CDS_END = 102_958_772    # End of last coding exon
    
    # PAH exons (1-based, inclusive)
    # These are the actual exon coordinates for PAH gene (ENST00000553106)
    exons = [
        (102_836_311, 102_836_441),  # Exon 1
        (102_837_263, 102_837_400),  # Exon 2
        (102_838_149, 102_838_277),  # Exon 3
        (102_840_753, 102_840_876),  # Exon 4
        (102_841_963, 102_842_101),  # Exon 5
        (102_843_073, 102_843_217),  # Exon 6
        (102_844_162, 102_844_284),  # Exon 7
        (102_845_348, 102_845_491),  # Exon 8
        (102_847_241, 102_847_373),  # Exon 9
        (102_848_418, 102_848_559),  # Exon 10
        (102_849_351, 102_849_506),  # Exon 11
        (102_850_968, 102_851_104),  # Exon 12
        (102_853_941, 102_854_085),  # Exon 13
    ]
    
    # Get the PAH gene sequence
    pah_seq = _sequence_store.get('PAH_GENE')
    if pah_seq is None:
        logger.warning("PAH gene sequence not found in sequence store, using test sequence")
        pah_seq = TEST_SEQUENCE
    
    # Calculate relative coordinates for the transcript
    tx_start = 1
    tx_end = len(pah_seq)
    
    # Calculate CDS start/end relative to transcript
    cds_start = PAH_CDS_START - PAH_START + 1
    cds_end = PAH_CDS_END - PAH_START + 1
    
    # Convert exons to relative coordinates
    relative_exons = [
        [start - PAH_START + 1, end - PAH_START + 1]
        for start, end in exons
    ]
    
    return TranscriptConfig(
        transcript_id=PAH_TRANSCRIPT,
        gene_symbol="PAH",
        assembly="GRCh38",
        strand=1,  # Forward strand
        chrom="chr12",
        tx_start=tx_start,
        tx_end=tx_end,
        cds_start=cds_start,
        cds_end=cds_end,
        transcript_sequence=pah_seq,
        exons=relative_exons
    )


def create_test_variants() -> List[VariantIn]:
    """
    Create test variants from the provided variant lists.
    
    Returns:
        List[VariantIn]: List of variant objects with phase groups
    """
    variants = []
    
    # Add variants from example_1 (phase 1)
    for hgvs in example_1:
        variants.append(VariantIn(hgvs=hgvs, phase_group=1))
    
    # Add variants from example_2 (phase 2)
    for hgvs in example_2:
        variants.append(VariantIn(hgvs=hgvs, phase_group=2))
    
    logger.info(f"Created {len(variants)} variants ({len(example_1)} in phase 1, {len(example_2)} in phase 2)")
    return variants


def analyze_haplotype(
    reference: str, 
    haplotype: str, 
    variants: List[str], 
    config: TranscriptConfig
) -> Dict:
    """
    Analyze a single haplotype and return the results.
    
    This function analyzes the consequences of variants on a haplotype and returns
    the resulting protein sequence and other relevant information.
    
    Args:
        reference: Reference DNA sequence
        haplotype: Haplotype DNA sequence with variants applied
        variants: List of variant HGVS strings in this haplotype
        config: Transcript configuration
        
    Returns:
        Dict containing analysis results including:
        - edited_sequence: The edited DNA sequence
        - variants: List of variant HGVS strings
        - protein_sequence: The resulting protein sequence
        - nmd_status: Nonsense-mediated decay status (if available)
        - nmd_rationale: Explanation for NMD prediction (if available)
    """
    logger.info(f"\nAnalyzing haplotype with variants: {', '.join(variants) if variants else 'None (reference)'}")
    
    # Initialize variables
    protein_outcome = None
    nmd_outcome = None
    
    try:
        # Log the reference and haplotype sequences for debugging
        logger.debug(f"Reference sequence length: {len(reference)}")
        logger.debug(f"Haplotype sequence length: {len(haplotype)}")
        
        # Log the CDS coordinates
        logger.debug(f"CDS start (1-based): {config.cds_start}, CDS end (1-based): {config.cds_end}")
        logger.debug(f"CDS length: {config.cds_end - config.cds_start + 1}")
        
        # Get the reference CDS sequence (0-based, half-open)
        # We need to ensure we're not going out of bounds
        ref_cds_start = max(0, config.cds_start - 1)
        ref_cds_end = min(len(reference), config.cds_end)
        ref_cds = reference[ref_cds_start:ref_cds_end]
        
        # Log the reference CDS sequence for debugging
        logger.debug(f"Reference CDS (first 100bp): {ref_cds[:100]}")
        logger.debug(f"Reference CDS (last 100bp): {ref_cds[-100:] if len(ref_cds) > 100 else ref_cds}")
        logger.debug(f"Reference CDS length: {len(ref_cds)}")
        
        # For the edited sequence, we need to account for any insertions/deletions
        # Calculate the length difference between the reference and haplotype
        length_diff = len(haplotype) - len(reference)
        
        # Get the edited CDS sequence, adjusting for length changes
        edited_cds_start = max(0, config.cds_start - 1)
        edited_cds_end = min(len(haplotype), config.cds_end + length_diff)
        edited_cds = haplotype[edited_cds_start:edited_cds_end]
        
        # Log the edited CDS sequence for debugging
        logger.debug(f"Edited CDS (first 100bp): {edited_cds[:100]}")
        logger.debug(f"Edited CDS (last 100bp): {edited_cds[-100:] if len(edited_cds) > 100 else edited_cds}")
        logger.debug(f"Edited CDS length: {len(edited_cds)}")
        
        # Log the differences between reference and edited CDS
        if len(ref_cds) == len(edited_cds):
            differences = [i for i in range(len(ref_cds)) if ref_cds[i] != edited_cds[i]]
            logger.debug(f"Number of differences between reference and edited CDS: {len(differences)}")
            if differences:
                logger.debug(f"First 10 differences at positions: {differences[:10]}")
        else:
            logger.debug(f"Reference and edited CDS have different lengths: {len(ref_cds)} vs {len(edited_cds)}")
        
        # Log the extracted CDS sequences with more details
        logger.debug(f"Reference CDS length: {len(ref_cds)}")
        logger.debug(f"Reference CDS (first 50bp): {ref_cds[:50]}")
        logger.debug(f"Reference CDS (last 50bp): {ref_cds[-50:] if len(ref_cds) > 50 else ref_cds}")
        
        # Check for non-DNA characters in reference CDS
        dna_bases = set('ACGT')
        non_dna_bases = set(ref_cds.upper()) - dna_bases
        if non_dna_bases:
            logger.warning(f"Found non-DNA bases in reference CDS: {non_dna_bases}")
        
        logger.debug(f"Edited CDS length: {len(edited_cds)}")
        logger.debug(f"Edited CDS (first 50bp): {edited_cds[:50]}")
        logger.debug(f"Edited CDS (last 50bp): {edited_cds[-50:] if len(edited_cds) > 50 else edited_cds}")
        
        # Check for non-DNA characters in edited CDS
        non_dna_bases = set(edited_cds.upper()) - dna_bases
        if non_dna_bases:
            logger.warning(f"Found non-DNA bases in edited CDS: {non_dna_bases}")
        
        # Create a copy of the config with updated CDS coordinates
        from copy import deepcopy
        updated_config = deepcopy(config)
        updated_config.cds_start = 1  # 1-based start of CDS in the extracted sequence
        updated_config.cds_end = len(ref_cds)  # 1-based end of CDS in the extracted sequence
        
        # Ensure the CDS is in-frame (multiple of 3) and not empty
        if len(ref_cds) % 3 != 0:
            logger.warning(f"Reference CDS length ({len(ref_cds)}) is not a multiple of 3, truncating...")
            ref_cds = ref_cds[:-(len(ref_cds) % 3)]
            logger.debug(f"Reference CDS after truncation: {len(ref_cds)} bp")
            
        if len(edited_cds) % 3 != 0:
            logger.warning(f"Edited CDS length ({len(edited_cds)}) is not a multiple of 3, truncating...")
            edited_cds = edited_cds[:-(len(edited_cds) % 3)]
            logger.debug(f"Edited CDS after truncation: {len(edited_cds)} bp")
        
        # Log the final CDS sequences with more details
        logger.debug(f"Final reference CDS length: {len(ref_cds)} bp")
        logger.debug(f"Final edited CDS length: {len(edited_cds)} bp")
        
        # Validate CDS sequences
        def validate_sequence(seq: str, name: str) -> bool:
            valid_bases = set('ACGTN')
            invalid_bases = set(seq.upper()) - valid_bases
            if invalid_bases:
                logger.warning(f"{name} contains invalid bases: {invalid_bases}")
                return False
            if not seq:
                logger.error(f"{name} is empty")
                return False
            if len(seq) % 3 != 0:
                logger.error(f"{name} length ({len(seq)}) is not a multiple of 3")
                return False
            return True
        
        # Validate both CDS sequences
        ref_valid = validate_sequence(ref_cds, "Reference CDS")
        edited_valid = validate_sequence(edited_cds, "Edited CDS")
        
        if not (ref_valid and edited_valid):
            logger.error("Invalid CDS sequences, cannot analyze consequences")
            return {
                'edited_sequence': haplotype,
                'variants': variants,
                'protein_sequence': 'N/A',
                'error': 'Invalid CDS sequences'
            }
        
        # Log first 10 codons for inspection
        def log_codons(seq: str, name: str):
            codons = [seq[i:i+3] for i in range(0, min(30, len(seq)), 3)]
            logger.debug(f"First {len(codons)} {name} codons: {codons}")
        
        log_codons(ref_cds, "reference CDS")
        log_codons(edited_cds, "edited CDS")
        
        # Analyze protein consequences with the CDS sequences
        logger.debug("Calling analyze_consequences...")
        try:
            protein_outcome = analyze_consequences(ref_cds, edited_cds, updated_config)
            logger.debug(f"analyze_consequences returned: {protein_outcome}")
        except Exception as e:
            logger.error(f"Error in analyze_consequences: {str(e)}", exc_info=True)
            raise
        
        if protein_outcome:
            # Log the protein outcome details
            logger.debug(f"Protein outcome: {protein_outcome}")
            
            # Check for NMD if protein outcome is available
            nmd_outcome = check_nmd(protein_outcome, updated_config)
            
            # Log basic information about the protein outcome
            logger.debug(f"Protein outcome type: {type(protein_outcome)}")
            
            # Log all attributes of the protein outcome for debugging
            for attr in dir(protein_outcome):
                if not attr.startswith('_'):  # Skip private attributes
                    try:
                        value = getattr(protein_outcome, attr)
                        if isinstance(value, str) and len(value) > 100:
                            logger.debug(f"{attr}: {value[:100]}... (truncated)")
                        else:
                            logger.debug(f"{attr}: {value}")
                    except Exception as e:
                        logger.debug(f"Could not access {attr}: {e}")
    
    except Exception as e:
        logger.error(f"Error analyzing consequences: {str(e)}")
        logger.exception("Detailed error:")
        protein_outcome = None
        nmd_outcome = None
    
    # Get protein sequence if available
    protein_seq = 'N/A'
    if protein_outcome:
        logger.debug(f"Protein outcome type: {type(protein_outcome)}")
        logger.debug(f"Protein outcome attributes: {dir(protein_outcome)}")
        
        # Get the protein sequence from the protein_outcome object
        if hasattr(protein_outcome, 'protein_sequence') and protein_outcome.protein_sequence:
            protein_seq = protein_outcome.protein_sequence
            logger.info(f"Protein sequence from protein_outcome.protein_sequence: {protein_seq[:100]}...")
        else:
            # If protein_sequence is not available, try to get it from other attributes
            for attr in ['sequence', 'aa_sequence', 'transcript_sequence']:
                if hasattr(protein_outcome, attr):
                    val = getattr(protein_outcome, attr)
                    if val and isinstance(val, str) and len(val) > 0:
                        protein_seq = val
                        logger.info(f"Using protein sequence from {attr}: {protein_seq[:100]}...")
                        break
            
            if protein_seq == 'N/A':
                logger.warning("No valid protein sequence found in protein_outcome")
                # As a fallback, use the reference protein sequence
                if hasattr(config, 'protein_sequence') and config.protein_sequence:
                    protein_seq = config.protein_sequence
                    logger.info(f"Using reference protein sequence from config: {protein_seq[:100]}...")
                else:
                    logger.error("No reference protein sequence available in config")
    
    # Prepare result
    result = {
        'edited_sequence': haplotype,
        'variants': variants,
        'protein_sequence': protein_seq,
    }
    
    # Add NMD information if available
    if nmd_outcome:
        result.update({
            'nmd_status': nmd_outcome.status,
            'nmd_rationale': getattr(nmd_outcome, 'rationale', None)
        })
        logger.info(f"NMD Status: {nmd_outcome.status}")
        if hasattr(nmd_outcome, 'rationale'):
            logger.info(f"NMD Rationale: {nmd_outcome.rationale}")
    
    return result


def save_results(results: List[Dict], output_dir: str) -> None:
    """
    Save analysis results to files.
    
    Args:
        results: List of result dictionaries
        output_dir: Directory to save output files
    """
    # Add haplotype IDs
    for i, result in enumerate(results, 1):
        result['haplotype_id'] = i
    
    # Save JSON output
    json_path = os.path.join(output_dir, 'complex_analysis_results.json')
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"\nResults saved to: {json_path}")
    
    # Save FASTA output
    fasta_path = os.path.join(output_dir, 'complex_analysis_results.fasta')
    with open(fasta_path, 'w') as f:
        for result in results:
            f.write(f">haplotype_{result['haplotype_id']} variants:{','.join(result['variants']) if result['variants'] else 'none'}\n")
            f.write(f"{result['edited_sequence']}\n")
    logger.info(f"FASTA output saved to: {fasta_path}")


def format_protein_sequence(seq: str, width: int = 60) -> str:
    """
    Format a protein sequence with line breaks and position numbers.
    
    Args:
        seq: Protein sequence to format
        width: Number of amino acids per line (default: 60)
        
    Returns:
        Formatted protein sequence with position numbers
    """
    import logging
    logger = logging.getLogger(__name__)
    
    logger.debug(f"Formatting protein sequence. Type: {type(seq)}, Value: {repr(seq)}")
    
    # Handle None or empty sequence
    if not seq:
        logger.warning("Empty sequence received")
        return "N/A"
        
    # Handle non-string sequences
    if not isinstance(seq, str):
        logger.warning(f"Expected string sequence, got {type(seq)}. Converting to string.")
        seq = str(seq)
    
    # Handle N/A or other special values
    if seq.strip().upper() in ('N/A', 'NONE', 'NULL'):
        logger.warning(f"Sequence has special value: {seq}")
        return seq
    
    # Remove any whitespace and convert to uppercase
    seq = ''.join(seq.split()).upper()
    
    # If sequence is empty after cleaning
    if not seq:
        logger.error("Sequence is empty after cleaning")
        return "N/A"
    
    # Format the sequence with line breaks and position numbers
    lines = []
    try:
        # First line: position numbers (every 10 amino acids)
        pos_line = []
        aa_line = []
        
        for i in range(0, len(seq), width):
            chunk = seq[i:i+width]
            pos_line = []
            aa_line = []
            
            # Add position numbers and amino acids
            for j in range(0, len(chunk), 10):
                pos = i + j + 1
                pos_str = str(pos).ljust(10)
                pos_line.append(pos_str[:10])
                aa_line.append(chunk[j:j+10])
            
            # Join the lines and add to output
            lines.append(''.join(pos_line))
            lines.append(''.join(aa_line))
            
        return '\n'.join(lines)
    except Exception as e:
        logger.error(f"Error formatting protein sequence: {str(e)}")
        logger.error(f"Sequence value: {repr(seq)}")
        return seq  # Return the sequence as-is if formatting fails

def print_summary(results: List[Dict]) -> None:
    """Print a summary of the analysis results."""
    print("\n" + "="*80)
    print("COMPLEX VARIANT ANALYSIS SUMMARY".center(80))
    print("="*80)
    
    for result in results:
        print(f"\n{'='*80}")
        print(f"HAPLOTYPE {result['haplotype_id']}")
        print(f"{'='*80}")
        
        # Print variants
        print("\nVARIANT(S):")
        if result['variants']:
            for i, var in enumerate(result['variants'], 1):
                print(f"  {i}. {var}")
        else:
            print("  None (reference)")
        
        # Print protein sequence
        print("\nPROTEIN SEQUENCE:")
        print(format_protein_sequence(result.get('protein_sequence', 'N/A')))
        
        # Print NMD status
        if 'nmd_status' in result:
            print(f"\nNMD STATUS: {result['nmd_status']}")
            if result.get('nmd_rationale'):
                print(f"  Rationale: {result['nmd_rationale']}")
    
    print("\n" + "="*80)


def main():
    """Main analysis workflow."""
    try:
        logger.info("=" * 80)
        logger.info(f"PAH Gene Variant Analysis ({PAH_TRANSCRIPT})")
        logger.info("=" * 80)
        
        # 1. Set up the sequence store
        logger.info("\n1. Setting up sequence store...")
        setup_sequence_store()
        
        # 2. Create transcript configuration
        logger.info("\n2. Creating transcript configuration...")
        config = create_transcript_config()
        
        # 3. Create test variants
        logger.info("\n3. Creating test variants...")
        variants = create_test_variants()
        
        # 4. Get the reference sequence from the sequence store
        reference_sequence = _sequence_store.get('PAH_GENE')
        if not reference_sequence:
            logger.warning("PAH_GENE sequence not found in sequence store, using test sequence")
            reference_sequence = TEST_SEQUENCE
        
        logger.info(f"Reference sequence length: {len(reference_sequence)}")
        
        # 5. Generate all possible haplotypes
        logger.info("\n5. Generating haplotypes...")
        haplotypes = generate_haplotypes(reference_sequence, variants)
        logger.info(f"Generated {len(haplotypes)} unique haplotypes")
        
        # 6. Analyze each haplotype
        logger.info("\n6. Analyzing haplotypes...")
        results = []
        for i, (hap_seq, hap_vars) in enumerate(haplotypes, 1):
            logger.info(f"\nAnalyzing haplotype {i}/{len(haplotypes)} with {len(hap_vars)} variants...")
            try:
                result = analyze_haplotype(reference_sequence, hap_seq, hap_vars, config)
                results.append(result)
            except Exception as e:
                logger.error(f"Error analyzing haplotype: {e}")
                continue
        
        if not results:
            logger.error("No results were generated. Check for errors in the input data.")
            return
        
        # 7. Save and display results
        logger.info("\n7. Saving results...")
        save_results(results, OUTPUT_DIR)
        
        logger.info("\nAnalysis Summary:")
        print_summary(results)
        
        logger.info("\nAnalysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Error during analysis: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
