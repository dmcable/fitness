from pyfaidx import Fasta
import pandas as pd

if(False):
    path = 'input_data/genome.fa'
    genome = Fasta(path)

    chrom = '1'
    pos = 100000

    # get a single base (on the + strand)
    genome[chrom][pos].seq

    # get a stretch of bases (end position is 0-based)
    genome[chrom][pos-1:pos+2].seq

def get_allele_mut_rates()
allele_mut_rates = pd.read_csv("input_data/rates.txt", delimiter=r"\s+")

def get_surrounding_genome(chrom, pos):
    return genome[chrom][pos-1:pos+2].seq

def get_muation_rate(chrom, pos, from_char, to_char):
    base_genome = get_surrounding_genome(chrom, pos)
    from_seq = base_genome[0] + from_char + base_genome[2]
    to_seq = base_genome[0] + to_char + base_genome[2]
    from_match = allele_mut_rates.loc[allele_mut_rates['from'] == from_seq]
    match = from_match.loc[from_match['to'] == to_seq]
    return match['mu_snp'].iloc[0]
