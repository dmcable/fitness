import pdb
from FitnessEstimator import FitnessEstimator
from load_data import load_data

if __name__ == '__main__':
    gene_df_ptv, gene_df_missense = load_data()
    estimator_ptv = FitnessEstimator(gene_df_ptv)
    #estimator_ptv.calculate_posterior_mean(gene_df_ptv)
