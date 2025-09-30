"""Script to demonstrate complex variant phasing with hgvs2seq."""
import sys
import os
import json
import logging
from typing import List, Dict, Tuple

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import hgvs2seq components
from hgvs2seq.data_provider import _sequence_store
from hgvs2seq.models import TranscriptConfig, VariantIn
from hgvs2seq.parse import parse_and_normalize_variants
from hgvs2seq.apply.batch import apply_variants_in_batch
from hgvs2seq.consequence.cds import analyze_consequences
from hgvs2seq.consequence.nmd import check_nmd

# Create a test transcript sequence with a valid CDS (must be multiple of 3)
# This is a made-up sequence with a start (ATG) and stop (TAA) codon
test_sequence = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTAATAA'

# Add the test sequence to the in-memory store
_sequence_store['TEST_TRANSCRIPT'] = test_sequence
_sequence_store['chr1'] = 'N' * 1000 + test_sequence + 'N' * 1000

# Print debug info
print(f"Added test sequence with length: {len(test_sequence)}")
print(f"Available sequence IDs: {list(_sequence_store.keys())}")

def generate_haplotypes(reference: str, variants: List[VariantIn]) -> List[Tuple[str, List[str]]]:
    """Generate all possible haplotypes from a list of variants with phase groups."""
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
    for phase, phase_vars in phase_groups.items():
        new_haplotypes = []
        for hap_seq, hap_vars in haplotypes:
            # Create two haplotypes for this phase (0 = reference, 1 = alternate)
            for allele in [0, 1]:
                new_hap = hap_seq.copy()
                new_vars = hap_vars.copy()
                
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
                        print(f"Error processing variant {var.hgvs}: {e}")
                
                new_haplotypes.append((new_hap, new_vars))
        
        # Update haplotypes for next phase
        haplotypes = new_haplotypes
    
    # Convert back to strings and remove duplicates
    unique_haplotypes = {}
    for hap_seq, hap_vars in haplotypes:
        hap_str = ''.join(hap_seq)
        if hap_str not in unique_haplotypes:
            unique_haplotypes[hap_str] = sorted(hap_vars)
    
    return [(hap, vars) for hap, vars in unique_haplotypes.items()]

def main():
    # 1. Create a test transcript configuration
    cds_start = 1
    cds_end = 138  # Just before the stop codon
    
    config = TranscriptConfig(
        transcript_id="TEST_TRANSCRIPT",
        gene_symbol="TEST_GENE",
        assembly="GRCh38",
        strand=1,
        chrom="chr1",
        tx_start=1,
        tx_end=len(test_sequence),
        cds_start=cds_start,
        cds_end=cds_end,
        transcript_sequence=test_sequence,
        exons=[[1, len(test_sequence)]]  # Single exon for simplicity
    )
    
    # 2. Create test variants with different phase groups
    # Format: (hgvs, phase_group, description)
    variant_data = [
        ("TEST_TRANSCRIPT:c.5G>A", 1, "Missense variant in phase 1"),
        ("TEST_TRANSCRIPT:c.10G>T", 1, "Missense variant in phase 1 (co-inherited with c.5G>A)"),
        ("TEST_TRANSCRIPT:c.15C>A", 2, "Missense variant in phase 2"),
        ("TEST_TRANSCRIPT:c.20T>G", 2, "Missense variant in phase 2 (co-inherited with c.15C>A)")
    ]
    
    variants = [VariantIn(hgvs=hgvs, phase_group=phase) 
               for hgvs, phase, _ in variant_data]
    
    try:
        # 3. Parse and normalize variants
        print("Parsing and normalizing variants...")
        norm_variants = parse_and_normalize_variants(variants, config)
        
        # 4. Generate all possible haplotypes
        print("\n=== Generating Haplotypes ===")
        haplotypes = generate_haplotypes(test_sequence, norm_variants)
        
        # 5. Process each haplotype
        print("\n=== Analyzing Haplotypes ===")
        results = []
        
        for hap_id, (hap_seq, hap_vars) in enumerate(haplotypes):
            print(f"\nHaplotype {hap_id + 1}:")
            print(f"Sequence: {hap_seq}")
            print(f"Variants: {', '.join(hap_vars) if hap_vars else 'None (reference)'}")
            
            # Analyze consequences
            protein_outcome = analyze_consequences(test_sequence, hap_seq, config)
            nmd_outcome = check_nmd(protein_outcome, config) if protein_outcome else None
            
            # Get protein sequence
            protein_seq = protein_outcome.protein_sequence if (protein_outcome and hasattr(protein_outcome, 'protein_sequence')) else None
            
            # Add to results
            results.append({
                'haplotype_id': hap_id + 1,
                'variants': hap_vars,
                'edited_sequence': hap_seq,
                'protein_sequence': protein_seq,
                'nmd_status': nmd_outcome.status if nmd_outcome else None,
                'nmd_rationale': getattr(nmd_outcome, 'rationale', None) if nmd_outcome else None
            })
            
            # Print summary
            print(f"Protein: {protein_seq}")
            if nmd_outcome:
                print(f"NMD Status: {nmd_outcome.status}")
                if hasattr(nmd_outcome, 'rationale'):
                    print(f"NMD Rationale: {nmd_outcome.rationale}")
        
        # 6. Save results
        output_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
        os.makedirs(output_dir, exist_ok=True)
        
        # Save JSON output
        json_path = os.path.join(output_dir, 'complex_test_results.json')
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to: {json_path}")
        
        # Save FASTA output
        fasta_path = os.path.join(output_dir, 'complex_test_results.fasta')
        with open(fasta_path, 'w') as f:
            for result in results:
                f.write(f">haplotype_{result['haplotype_id']} variants:{','.join(result['variants']) if result['variants'] else 'none'}\n")
                f.write(f"{result['edited_sequence']}\n")
        print(f"FASTA output saved to: {fasta_path}")
        
        # Print a summary
        print("\n=== Final Summary ===")
        for result in results:
            print(f"\nHaplotype {result['haplotype_id']}:")
            print(f"Variants: {', '.join(result['variants']) if result['variants'] else 'None (reference)'}")
            print(f"Protein: {result['protein_sequence']}")
            if result['nmd_status']:
                print(f"NMD: {result['nmd_status']} - {result.get('nmd_rationale', '')}")
        
        # Print variant descriptions
        print("\n=== Variant Descriptions ===")
        for hgvs, _, desc in variant_data:
            print(f"{hgvs}: {desc}")
        
    except Exception as e:
        logger.error(f"Error running hgvs2seq analysis: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
