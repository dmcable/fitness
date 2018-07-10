import pandas as pd
import pdb
import numpy as np
import pickle


# creates a string lookup key by concatenating the position, reference, and alt allele
def get_lookup_key(pos, ref, alt):
    return pos + ref + alt


# returns variant_scores, a mapping from lookup_key to primate_score quartile per gene
def get_variant_quartiles():
    primate_path = 'input_data/primate_scores.txt'
    variant_scores = pd.DataFrame(columns=['gene', 'score'])
    first_line = True
    gene_dict = {}
    count = 0
    with open(primate_path, 'r') as f:
        for line in f:
            count += 1
            if(count % 10000 == 0):
                print(count)
            if(count > 100000):
                break
            if(not first_line):
                elems = line.split('\t')
                lookup_key = get_lookup_key(elems[1], elems[2], elems[3])
                curr_score = float(elems[5][:-1])
                new_pt = pd.DataFrame(data = [{'gene':elems[4], 'score':curr_score}], index = [lookup_key])
                variant_scores = variant_scores.append(new_pt)
                if(elems[4] not in gene_dict):
                    gene_dict[elems[4]] = []
                gene_dict[elems[4]].append(curr_score)
            else:
                first_line = False
    gene_quartile_fun = {}
    for gene, value_list in gene_dict.iteritems():
        drop_gene = len(value_list) < 10
        gene_values = np.array(value_list)
        cutoffs = [np.percentile(gene_values, 25), np.percentile(gene_values, 50), np.percentile(gene_values, 75)]

        def get_quartile(value):
            if(drop_gene):
                return -1
            elif(value < cutoffs[0]):
                return '0'
            elif(value < cutoffs[1]):
                return '1'
            elif(value < cutoffs[2]):
                return '2'
            return '3'

        gene_quartile_fun[gene] = get_quartile
    variant_quartiles = {}
    for key in variant_scores.index:
        gene = variant_scores.loc[key, 'gene']
        score = variant_scores.loc[key, 'score']
        quartile = gene_quartile_fun[gene](score)
        variant_quartiles[key] = quartile
    with open('saved_data/variant_quartiles.pickle', 'wb') as handle:
        pickle.dump(variant_quartiles, handle, protocol=pickle.HIGHEST_PROTOCOL)
    variant_scores.to_pickle('saved_data/variant_scores.pkl')
    return variant_scores
    # Load data (deserialize)
    #with open('filename.pickle', 'rb') as handle:
        #unserialized_data = pickle.load(handle)
    pdb.set_trace()
