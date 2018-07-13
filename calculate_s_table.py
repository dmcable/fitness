import pdb
from FitnessEstimator import FitnessEstimator
from load_data import load_data
from load_data import load_syn
import pandas as pd


if __name__ == '__main__':
    type = 'syn'
    assert type in ['syn', 'mis', 'ptv']
    if(type == 'mis'):
        _, gene_df = load_data()
    elif(type == 'ptv'):
        gene_df, _ = load_data()
    else:
        gene_df = load_syn()
    estimator = FitnessEstimator(gene_df)
    estimator.calculate_posterior_mean(gene_df)
    gene_df.to_pickle('saved_data/gene_df_' + type + '_means.pkl')
