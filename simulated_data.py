import numpy as np
import pdb
import pandas as pd


def simulate_data(gene_df, estimator):
    N_genomes = 60706
    simulated_df = pd.DataFrame.copy(gene_df)
    simulated_df['s_het'] = estimator.sample_prior(len(simulated_df))
    n_mean = simulated_df['total_mutations'] / simulated_df['s_het']
    simulated_df['allele_count'] = np.random.poisson(n_mean)
    gene_filter = simulated_df['allele_count'] <= 0.001 * N_genomes
    simulated_df = simulated_df[gene_filter]
    return simulated_df
