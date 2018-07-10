import pdb
from FitnessEstimator import FitnessEstimator
from load_data import load_data
from load_data import load_syn
from simulated_data import simulate_data
import pandas as pd
import numpy as np


def evaluate_prior_method(gene_df):
    estimator = FitnessEstimator(gene_df, run_stuff = True)
    simulated_df = simulate_data(gene_df.sample(1000), estimator)
    estimator = FitnessEstimator(simulated_df, run_stuff=True)
    N_genes = [1000] #1000, 3000
    #trials = [2, 2, 2, 2, 2, 2]
    #trials = [50, 50, 50, 20, 10, 5]
    trials = [10]
    methods = ['iterative']
    results = np.zeros((len(N_genes), max(trials), len(methods), 2))
    for t in range(max(trials)):
        print(t)
        for method_index in range(len(methods)):
            method = methods[method_index]
            for cond in range(len(N_genes)):
                print(cond)
                if(t < trials[cond]):
                    N_gene = N_genes[cond]
                    simulated_df = simulate_data(gene_df.sample(N_gene), estimator)
                    alpha, beta = estimator.estimate_prior(simulated_df, method=method, grid_size=10)
                    print(str(alpha) + " " + str(beta))
                    results[cond, t, method_index, 0] = alpha
                    results[cond, t, method_index, 1] = beta
    means = np.mean(results, axis=1)
    stds = np.std(results, axis=1)
    np.save('results/means2_' + methods[0] + '.npy', means)
    np.save('results/stds2_' + methods[0] + '.npy', stds)
    pdb.set_trace()


if __name__ == '__main__':
    #gene_df, gene_df_missense = load_data()
    gene_df = pd.read_pickle('saved_data/gene_df_strat.pkl')
    #simulated_df = pd.read_pickle('saved_data/simulated_data.pkl')
    #evaluate_prior_method(gene_df_ptv)
    #gene_df_ptv = gene_df_ptv.sample(1000)
    estimator = FitnessEstimator(gene_df)
    estimator.calculate_posterior_mean_strat(gene_df)
    #gene_df = gene_df[gene_df['s_mean'] > 0]
    #gene_df = gene_df[gene_df['allele_count'] > 0]
    gene_df.to_pickle('saved_data/gene_df_strat_means.pkl')
    pdb.set_trace()
    #simulated_df = simulate_data(gene_df_ptv, estimator)
    #estimator.estimate_prior(simulated_df, method='empirical', grid_size=10)
