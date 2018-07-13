import pandas as pd
import pdb
import numpy as np


# creates a string lookup key by concatenating the position, reference, and alt allele
def get_lookup_key(pos, ref, alt):
    return pos + ref + alt


# variant_scores is a dataframe mapping variant keys to gene and primate scores
# gene_scores is dictionary mapping gene name to a list of primate scores across its variants
def get_variant_quartiles():
    primate_path = 'input_data/primate_scores.txt'
    quartiles_list = []
    index_list = []
    with open(primate_path, 'r') as f:
        curr_gene = ""
        variant_scores = []
        variant_index_list = []
        for line in f:
            elems = line.split('\t')
            if(elems[0] != 'chrom'):
                pos, ref, alt, gene = (elems[1], elems[2], elems[3], elems[4])
                score = float(elems[5][:-1])
                lookup_key = get_lookup_key(pos, ref, alt)
                if(gene != curr_gene):
                    if(curr_gene != ''):
                        gene_quartile_fun = get_gene_quartile_fun(variant_scores)
                        index_list += variant_index_list
                        quartiles_list += map(gene_quartile_fun, variant_scores)
                    curr_gene = gene
                    variant_scores = []
                    variant_index_list = []
                variant_scores.append(score)
                variant_index_list.append(lookup_key)
    return pd.DataFrame({'quartile': quartiles_list}, index=index_list)


def get_gene_quartile_fun(variant_scores):
    gene_values = np.array(variant_scores)
    cutoffs = [
        np.percentile(gene_values, 25), np.percentile(gene_values, 50),
        np.percentile(gene_values, 75)
    ]

    def get_quartile(value):
        if(value < cutoffs[0]):
            return '0'
        elif(value < cutoffs[1]):
            return '1'
        elif(value < cutoffs[2]):
            return '2'
        return '3'

    return get_quartile


# returns variant_quartiles, a mapping from lookup_key to primate_score quartile per gene
def load_variant_quartiles():
    quartiles_file = 'saved_data/variant_quartiles_df.pkl'
    loaded_variant_quartiles = True
    if(loaded_variant_quartiles):
        variant_quartiles = pd.read_pickle(quartiles_file)
    else:
        variant_quartiles = get_variant_quartiles()
        variant_quartiles.to_pickle(quartiles_file)
    return variant_quartiles


if __name__ == '__main__':
    load_variant_quartiles()
    pdb.set_trace()
