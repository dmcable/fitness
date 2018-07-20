import pandas as pd
import pdb
from FitnessEstimator import FitnessEstimator
from VariantQuartiles import load_variant_quartiles
from load_data import create_gene_df_stratified
from load_data import get_coverage_filter


def load_df(strat_df_file):
    mut_rates_path = 'input_data/mutation_probabilities.xls'
    mutation_rates = pd.read_excel(mut_rates_path, index_col=1, sheet_name=1)
    allele_df = pd.read_pickle('saved_data/allele_df.pkl')
    strat_names = ['0', '1', '2', '3']
    variant_quartiles = load_variant_quartiles()
    gene_df = create_gene_df_stratified(
        allele_df, mutation_rates, strat_names, variant_quartiles
    )
    gene_df.to_pickle(strat_df_file)
    return gene_df


def load_means(strat_df_file, strat_means_file, gene_df=None):
    if(gene_df is None):
        gene_df = pd.read_pickle(strat_df_file)
    estimator = FitnessEstimator(gene_df)
    estimator.calculate_posterior_mean_strat(gene_df)
    gene_df.to_pickle(strat_means_file)
    return gene_df


def augment_means(strat_means_file, gene_df=None):
    if(gene_df is None):
        gene_df = pd.read_pickle(strat_means_file)
    gene_df_ptv = pd.read_pickle('saved_data/gene_df_ptv_means.pkl')
    gene_df['s_ptv'] = gene_df_ptv['s_mean']
    gene_df['ac_ptv'] = gene_df_ptv['allele_count']
    gene_df_syn = pd.read_pickle('saved_data/gene_df_syn_means.pkl')
    gene_df['s_syn'] = gene_df_syn['s_mean']
    gene_df['ac_syn'] = gene_df_syn['allele_count']
    total_missense = gene_df['ac_0'] + gene_df['ac_1'] + gene_df['ac_2'] + gene_df['ac_3']
    gene_df = gene_df[total_missense > 1]
    gene_filter = get_coverage_filter(gene_df.index)
    gene_df = gene_df[gene_filter]
    gene_df.to_pickle('saved_data/gene_df_strat_augmented.pkl')
    return gene_df


def load_strat():
    strat_df_file = 'saved_data/gene_df_strat.pkl'
    strat_means_file = 'saved_data/gene_df_strat_means.pkl'
    modes = {'load_df', 'load_means', 'augment'}
    gene_df = None
    if('load_df' in modes):
        gene_df = load_df(strat_df_file)
    if('load_means' in modes):
        gene_df = load_means(strat_df_file, strat_means_file, gene_df)
    if('augment' in modes):
        gene_df = augment_means(strat_means_file, gene_df)
    return gene_df


# assumes modes is a continuous block
if __name__ == '__main__':
    gene_df = load_strat()
    pdb.set_trace()
