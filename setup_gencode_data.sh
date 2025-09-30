#!/bin/bash

# Create data directory if it doesn't exist
mkdir -p data

# Change to data directory
cd data

# Download chromosome 12 from GENCODE (GRCh38)
echo "Downloading chromosome 12 FASTA..."
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz

# Extract the FASTA
echo "Extracting chromosome 12..."
gunzip -f chr12.fa.gz
mv chr12.fa GRCh38.chr12.fasta

# Download comprehensive gene annotation
echo "Downloading comprehensive annotation..."
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz

# For now, create a minimal chr12 annotation file since the full extraction is problematic
zgrep '^chr12' data/gencode.v42.annotation.gtf.gz > data/GenCode42.chr12.gtf