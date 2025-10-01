import os
import unittest
from hgvs2seq.pyhgvs.utils import read_transcripts_refgene

class TestRefGeneReader(unittest.TestCase):
    def setUp(self):
        # Create a temporary refGene file for testing
        self.test_refgene = "test_refgene.txt"
        with open(self.test_refgene, 'w') as f:
            f.write("""585	NM_001346897.1	chr1	+	67091138	67216822	67091138	67216822	3	67091138,67105453,67216428,	67091551,67105573,67216822,	0	CDK11B	cmpl	cmpl	0,0,0,
586	NM_001346897.1	chr1	+	67091138	67216822	67091138	67216822	3	67091138,67105453,67216428,	67091551,67105573,67216822,	0	CDK11B	cmpl	cmpl	0,0,0,""")
    
    def tearDown(self):
        # Clean up the test file
        if os.path.exists(self.test_refgene):
            os.remove(self.test_refgene)
    
    def test_read_transcripts_refgene(self):
        with open(self.test_refgene) as f:
            transcripts = read_transcripts_refgene(f)
            
        # Should find at least one transcript (duplicates are removed by set)
        self.assertGreaterEqual(len(transcripts), 1)
        
        # Get the first transcript
        tx_id = next(iter(transcripts))
        tx = transcripts[tx_id]
        
        # Check basic properties
        self.assertEqual(tx.name, "NM_001346897")
        self.assertEqual(tx.version, 1)
        self.assertEqual(tx.gene.name, "CDK11B")
        
        # Check position
        self.assertEqual(tx.tx_position.chrom, "chr1")
        self.assertEqual(tx.tx_position.is_forward_strand, True)
        self.assertEqual(tx.tx_position.chrom_start, 67091138)  # 1-based start
        self.assertEqual(tx.tx_position.chrom_stop, 67216822)   # 1-based, end-inclusive
        
        # Check exons
        self.assertEqual(len(tx.exons), 3)
        
        # Check first exon
        exon1 = tx.exons[0]
        self.assertEqual(exon1.exon_number, 1)
        self.assertEqual(exon1.tx_position.chrom_start, 67091138)  # 1-based start
        self.assertEqual(exon1.tx_position.chrom_stop, 67091551)   # 1-based, end-inclusive

if __name__ == '__main__':
    unittest.main()
