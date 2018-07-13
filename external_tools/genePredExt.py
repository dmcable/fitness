
from collections import namedtuple

Gencode = namedtuple('Gencode', ['id', 'chrom', 'strand', 'tx_start',
    'tx_end', 'cds_start', 'cds_end', 'exons_n', 'exons',
    'score', 'name2', 'cds_start_stat', 'cds_end_stat', 'frames'])

def parse_genePredExt(line):
    ''' parse a genePredExt formatted line
    
    Parse a gencode line, which has been formatted in the genePredExt format,
    see https://genome.ucsc.edu/FAQ/FAQformat.html#format9.
    
    Args:
        line: genePredExt formatted line
    
    Yields:
        for each exon, yields chromosome and tuple of (start, end) positions
    '''
    tx_id, chrom, strand, tx_start, tx_end, cds_start, cds_end, exons_n, \
        exon_starts, exon_ends, score, name2, cds_start_stat, cds_end_stat, \
        frames = line.strip().split('\t')
    chrom = chrom.strip('chr').replace('M', 'MT')
    
    exon_starts = map(int, exon_starts.strip(',').split(','))
    exon_ends = map(int, exon_ends.strip(',').split(','))
    
    return Gencode(tx_id, chrom, strand, int(tx_start), int(tx_end),
        int(cds_start), int(cds_end), exons_n, list(zip(exon_starts, exon_ends)),
        score, name2, cds_start_stat, cds_end_stat, frames)
