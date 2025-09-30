
import pyhgvs
from pyfaidx import Fasta
from hgvs2seq.pyhgvs.utils import read_transcripts_gtf

# Read genome sequence using pyfaidx.
genome = Fasta('data/GRCh38.chr12.fa')

# Read RefSeq transcripts into a python dict.
with open('data/GenCode42.chr12.gtf') as infile:
    transcripts = read_transcripts_gtf(infile)

## Here is where the bug is. We are not getting a list of transcripts, which is a problem. 
## This comes from the change in the read_trascripts function in pyhgvs. Because the file format of gencode.v42.chr12.annotation.gtf
## is not the same as the one expected by pyhgvs, we get an empty list.

# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)

# Parse the HGVS name into genomic coordinates and alleles.
chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
    'ENST00000553106.6:c.1340C>A', genome, get_transcript=get_transcript)


print(chrom, offset, ref, alt)

# Format an HGVS name.
transcript = get_transcript('ENST00000553106.6')
hgvs_name = pyhgvs.format_hgvs_name(
    chrom, offset, ref, alt, genome, transcript)

print(hgvs_name)
