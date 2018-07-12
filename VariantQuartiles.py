import pandas as pd
import pdb
import numpy as np
import pickle


# creates a string lookup key by concatenating the position, reference, and alt allele
def get_lookup_key(pos, ref, alt):
    return pos + ref + alt


# variant_scores is a dataframe mapping variant keys to gene and primate scores
# gene_scores is dictionary mapping gene name to a list of primate scores across its variants
#@profile
def get_variant_quartiles():
    primate_path = 'input_data/primate_scores.txt'
    variant_quartiles = {}
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
                        for i in range(len(variant_index_list)):
                            variant_quartiles[variant_index_list[i]] = (
                                gene_quartile_fun(variant_scores[i])
                            )
                    curr_gene = gene
                    variant_scores = []
                    variant_index_list = []
                variant_scores.append(score)
                variant_index_list.append(lookup_key)
    return variant_quartiles


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
    quartiles_file = 'saved_data/variant_quartiles.pkl'
    loaded_variant_quartiles = False
    if(loaded_variant_quartiles):
        with open(quartiles_file, 'rb') as handle:
            variant_quartiles = pickle.load(handle)
    else:
        variant_quartiles = get_variant_quartiles()
        with open(quartiles_file, 'wb') as handle:
            pickle.dump(variant_quartiles, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return variant_quartiles


if __name__ == '__main__':
    load_variant_quartiles()
