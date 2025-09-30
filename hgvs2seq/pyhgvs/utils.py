"""
Helper functions.
"""

from __future__ import absolute_import
from __future__ import unicode_literals

import re
from collections import defaultdict

from .models import Exon
from .models import Position
from .models import Transcript


def read_refgene(infile):
    """
    Iterate through a refGene file.

    GenePred extension format:
    http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt

    Column definitions:
    0. uint undocumented id
    1. string name;             "Name of gene (usually transcript_id from GTF)"
    2. string chrom;                "Chromosome name"
    3. char[1] strand;              "+ or - for strand"
    4. uint txStart;                "Transcription start position"
    5. uint txEnd;                  "Transcription end position"
    6. uint cdsStart;               "Coding region start"
    7. uint cdsEnd;                 "Coding region end"
    8. uint exonCount;              "Number of exons"
    9. uint[exonCount] exonStarts;  "Exon start positions"
    10. uint[exonCount] exonEnds;   "Exon end positions"
    11. uint id;                    "Unique identifier"
    12. string name2;               "Alternate name (e.g. gene_id from GTF)"
    13. string cdsStartStat;        "enum('none','unk','incmpl','cmpl')"
    14. string cdsEndStat;          "enum('none','unk','incmpl','cmpl')"
    15. lstring exonFrames;         "Exon frame offsets {0,1,2}"
    """
    for line in infile:
        # Skip comments.
        if line.startswith('#'):
            continue
        row = line.rstrip('\n').split('\t')
        if len(row) != 16:
            raise ValueError(
                'File has incorrect number of columns '
                'in at least one line.')

        # Skip trailing ,
        exon_starts = list(map(int, row[9].split(',')[:-1]))
        exon_ends = list(map(int, row[10].split(',')[:-1]))
        exon_frames = list(map(int, row[15].split(',')[:-1]))
        exons = list(zip(exon_starts, exon_ends))

        yield {
            'chrom': row[2],
            'start': int(row[4]),
            'end': int(row[5]),
            'id': row[1],
            'strand': row[3],
            'cds_start': int(row[6]),
            'cds_end': int(row[7]),
            'gene_name': row[12],
            'exons': exons,
            'exon_frames': exon_frames
        }


def make_transcript(transcript_json):
    """
    Make a Transcript form a JSON object.
    """

    transcript_name = transcript_json['id']
    if '.' in transcript_name:
        name, version = transcript_name.split('.')
    else:
        name, version = transcript_name, None

    transcript = Transcript(
        name=name,
        version=int(version) if version is not None else None,
        gene=transcript_json['gene_name'],
        tx_position=Position(
            transcript_json['chrom'],
            transcript_json['start'],
            transcript_json['end'],
            transcript_json['strand'] == '+'),
        cds_position=Position(
            transcript_json['chrom'],
            transcript_json['cds_start'],
            transcript_json['cds_end'],
            transcript_json['strand'] == '+'))

    exons = transcript_json['exons']
    if not transcript.tx_position.is_forward_strand:
        exons = reversed(exons)

    for exon_number, (exon_start, exon_end) in enumerate(exons, 1):
        transcript.exons.append(
            Exon(transcript=transcript,
                 tx_position=Position(
                     transcript_json['chrom'],
                     exon_start,
                     exon_end,
                     transcript_json['strand'] == '+'),
                 exon_number=exon_number))

    return transcript

def read_transcripts(refgene_file):
    """
    Read all transcripts in a RefGene file.
    
    Args:
        refgene_file: A file-like object containing refGene data
        
    Returns:
        set: A set of Transcript objects
    """
    transcripts = set()
    for transcript_data in read_refgene(refgene_file):
        # Create a Transcript object from the dictionary
        transcript = make_transcript(transcript_data)
        transcripts.add(transcript)
    return transcripts


def _parse_gtf_attributes(attribute_string):
    """Parse the attributes string from a GTF file into a dictionary."""
    attributes = {}
    # Handle quoted values with semicolons inside them
    for match in re.finditer(r'([\w_]+)\s+"([^"]+)"', attribute_string):
        key, value = match.groups()
        attributes[key] = value
    return attributes


def read_transcripts_gtf(gtf_file):
    """
    Read transcripts from a GTF file.
    
    Args:
        gtf_file: A file-like object containing GTF data
        
    Returns:
        set: A set of Transcript objects
    """
    # First pass: collect all exons for each transcript
    transcript_exons = defaultdict(list)
    transcript_info = {}  # Store transcript metadata
    
    for line in gtf_file:
        # Skip comments and empty lines
        if line.startswith('#') or not line.strip():
            continue
            
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
            
        chrom, source, feature_type, start, end, score, strand, phase, attributes_str = fields[:9]
        start = int(start) - 1  # Convert to 0-based, half-open
        end = int(end)  # GTF is 1-based, end-inclusive; we'll keep as end-exclusive
        
        # Parse attributes
        attributes = _parse_gtf_attributes(attributes_str)
        
        transcript_id = attributes.get('transcript_id')
        if not transcript_id:
            continue
            
        gene_name = attributes.get('gene_name', transcript_id.split('.')[0])
        
        # Store transcript info if this is a transcript feature
        if feature_type == 'transcript':
            transcript_info[transcript_id] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand == '+',
                'gene_name': gene_name,
                'transcript_type': attributes.get('transcript_type', ''),
                'transcript_name': attributes.get('transcript_name', transcript_id),
            }
        # Store exon information
        elif feature_type == 'exon':
            exon_number = int(attributes.get('exon_number', 0))
            transcript_exons[transcript_id].append({
                'start': start,
                'end': end,
                'exon_number': exon_number,
                'strand': strand == '+',
            })
    
    # Second pass: create Transcript objects
    transcripts = set()
    
    for tx_id, info in transcript_info.items():
        # Skip if we don't have any exons for this transcript
        if tx_id not in transcript_exons:
            continue
        
        # Sort exons by position (5' to 3')
        exons = sorted(transcript_exons[tx_id], 
                      key=lambda e: e['start'], 
                      reverse=not info['strand'])
        
        # Create transcript position
        tx_position = Position(
            chrom=info['chrom'],
            chrom_start=info['start'],
            chrom_stop=info['end'],
            is_forward_strand=info['strand']
        )
        
        # For GTF, we'll use the transcript bounds as CDS for now
        # In a real implementation, you might want to parse CDS features
        cds_position = Position(
            chrom=info['chrom'],
            chrom_start=info['start'],
            chrom_stop=info['end'],
            is_forward_strand=info['strand']
        )
        
        # Create transcript
        tx_name, _, version = tx_id.partition('.')
        version = int(version) if version else None
        
        transcript = Transcript(
            name=tx_name,
            version=version,
            gene=info['gene_name'],
            tx_position=tx_position,
            cds_position=cds_position,
            is_default=('tag' in info and 'basic' in info['tag']),
            exons=[]
        )
        
        # Add exons to transcript
        for i, exon_info in enumerate(exons, 1):
            exon_position = Position(
                chrom=info['chrom'],
                chrom_start=exon_info['start'],
                chrom_stop=exon_info['end'],
                is_forward_strand=info['strand']
            )
            exon = Exon(
                transcript=transcript,
                tx_position=exon_position,
                exon_number=i
            )
            transcript.exons.append(exon)
            
        transcripts.add(transcript)
    
    return transcripts
