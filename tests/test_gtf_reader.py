import os
import unittest
from hgvs2seq.pyhgvs.utils import read_transcripts_gtf

class TestGTFReader(unittest.TestCase):    
    def setUp(self):
        # Create a temporary GTF file for testing
        self.test_gtf = "test.gtf"
        with open(self.test_gtf, 'w') as f:
            f.write("""# Test GTF file
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; tag "Ensembl_canonical";
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; tag "Ensembl_canonical";
chr1\tHAVANA\texon\t12613\t12721\t.\t+\t.\tgene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; tag "Ensembl_canonical";
""")
    
    def tearDown(self):
        # Clean up the test file
        if os.path.exists(self.test_gtf):
            os.remove(self.test_gtf)
    
    def test_read_transcripts_gtf(self):
        with open(self.test_gtf) as f:
            transcripts = read_transcripts_gtf(f)
            
        # Should find at least one transcript (duplicates are removed by set)
        self.assertGreaterEqual(len(transcripts), 1)
        
        # Get the first transcript
        tx = next(iter(transcripts))
        
        # Check basic properties
        self.assertEqual(tx.name, "ENST00000456328")
        self.assertEqual(tx.version, 2)
        self.assertEqual(tx.gene.name, "DDX11L1")
        self.assertTrue(tx.tx_position.is_forward_strand)
        
        # Check position
        self.assertEqual(tx.tx_position.chrom, "chr1")
        self.assertEqual(tx.tx_position.chrom_start, 11868)  # 0-based
        self.assertEqual(tx.tx_position.chrom_stop, 14409)   # 1-based, end-inclusive
        
        # Check exons
        self.assertEqual(len(tx.exons), 2)
        
        # Check first exon
        exon1 = tx.exons[0]
        self.assertEqual(exon1.exon_number, 1)
        self.assertEqual(exon1.tx_position.chrom_start, 11868)  # 0-based
        self.assertEqual(exon1.tx_position.chrom_stop, 12227)   # 1-based, end-inclusive
        
        # Check second exon
        exon2 = tx.exons[1]
        self.assertEqual(exon2.exon_number, 2)
        self.assertEqual(exon2.tx_position.chrom_start, 12612)  # 0-based
        self.assertEqual(exon2.tx_position.chrom_stop, 12721)   # 1-based, end-inclusive

if __name__ == '__main__':
    unittest.main()
