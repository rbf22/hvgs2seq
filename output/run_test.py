"""Script to set up test data and run hgvs2seq analysis directly."""
import sys
import os
import json
import logging
from typing import List

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
from hgvs2seq.io.jsonio import generate_json_output
from hgvs2seq.io.fasta import generate_fasta_output

# Create a test transcript sequence with a valid CDS (must be multiple of 3)
# This is a made-up sequence with a start (ATG) and stop (TAA) codon
test_sequence = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTAATAA'

# Add the test sequence to the in-memory store
_sequence_store['TEST_TRANSCRIPT'] = test_sequence
_sequence_store['chr1'] = 'N' * 1000 + test_sequence + 'N' * 1000

# Print debug info
print(f"Added test sequence with length: {len(test_sequence)}")
print(f"Available sequence IDs: {list(_sequence_store.keys())}")

def main():
    # 1. Create a test transcript configuration
    # Define CDS coordinates that match our test sequence
    # CDS starts at position 1 and ends before the stop codon (TAA at position 139-141)
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
    
    # 2. Create a test variant
    variant = VariantIn(hgvs="TEST_TRANSCRIPT:c.5G>A", phase_group=None)
    
    try:
        # 3. Parse and normalize the variant
        print("Parsing and normalizing variant...")
        norm_variants = parse_and_normalize_variants([variant], config)
        
        # 4. Apply the variant to the reference sequence
        print("Applying variant to reference sequence...")
        edited_sequences = apply_variants_in_batch(test_sequence, norm_variants)
        
        # 5. Analyze consequences for each haplotype
        print("Analyzing consequences...")
        results = []
        for hap_id, edited_seq in edited_sequences.items():
            protein_outcome = analyze_consequences(test_sequence, edited_seq, config)
            nmd_outcome = check_nmd(protein_outcome, config)
            
            results.append({
                'haplotype_id': hap_id,
                'edited_sequence': edited_seq,
                'protein_sequence': protein_outcome.protein_sequence if (protein_outcome and protein_outcome.protein_sequence) else "No protein sequence generated",
                'nmd_status': nmd_outcome.status if nmd_outcome else None,
                'nmd_rationale': getattr(nmd_outcome, 'rationale', 'No rationale provided') if nmd_outcome else None
            })
        
        # 6. Save results
        output_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
        os.makedirs(output_dir, exist_ok=True)
        
        # Save JSON output
        json_path = os.path.join(output_dir, 'test_results.json')
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to: {json_path}")
        
        # Save FASTA output
        fasta_path = os.path.join(output_dir, 'test_results.fasta')
        with open(fasta_path, 'w') as f:
            for result in results:
                f.write(f">haplotype_{result['haplotype_id']}\n{result['edited_sequence']}\n")
        print(f"FASTA output saved to: {fasta_path}")
        
        # Print a summary
        print("\n=== Summary ===")
        for result in results:
            print(f"\nHaplotype {result['haplotype_id']}:")
            print(f"Edited sequence: {result['edited_sequence']}")
            print(f"Protein sequence: {result['protein_sequence']}")
            if 'nmd_status' in result:
                print(f"NMD status: {result['nmd_status']}")
                if 'nmd_rationale' in result:
                    print(f"NMD rationale: {result['nmd_rationale']}")
            
    except Exception as e:
        logger.error(f"Error running hgvs2seq analysis: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
