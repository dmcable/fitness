import pdb
from FitnessEstimator import FitnessEstimator
from load_data import load_data
from load_data import load_syn
from load_data import load_global_allele_df
from analyze_strat import load_strat
import pandas as pd


def load_gene_df(type):
    type = 'syn'
    gene_df_file = 'saved_data/gene_df_' + type + '_means.pkl'
    pre_loaded = ['']
    assert type in ['syn', 'mis', 'ptv']
    if type in pre_loaded:
        gene_df = pd.read_pickle(gene_df_file)
    else:
        if(type == 'mis'):
            _, gene_df = load_data()
        elif(type == 'ptv'):
            gene_df, _ = load_data()
        elif(type == 'syn'):
            gene_df = load_syn()
        estimator = FitnessEstimator(gene_df)
        estimator.calculate_posterior_mean(gene_df)
        gene_df.to_pickle(gene_df_file)
    return gene_df


# performs our main analysis for the steady state method: we load in the gene_df
# for synonomys and ptv, and then we load in the stratified gene_df from missense
# primate scores. Then we calculate shet scores for all of these conditions.
if __name__ == '__main__':
    load_global_allele_df()
    types = ['syn', 'mis', 'ptv']
    for type in types:
        print('loading gene_df of type ' + type)
        load_gene_df(type)
    gene_df = load_strat()
