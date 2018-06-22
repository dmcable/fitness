from scipy.special import kv as bessel
import numpy as np
import matplotlib.pyplot as plt
import pdb
import math
import scipy.stats
import plotting


def compute_log_bessel(n, x):
    bessel_val = bessel(n, x)
    if(bessel_val == 0.0):
        return -np.inf
    log_bessel_val = np.log(bessel_val)
    if(np.isnan(log_bessel_val) or np.isinf(log_bessel_val)):
        raise ValueError
    return log_bessel_val


def assert_floats(float_list):
    for elem in float_list:
        assert(isinstance(elem, float))


class FitnessEstimator(object):

    @staticmethod
    def log_likelihood_analytical(alpha, beta, allele_count, total_mutations):
        assert_floats([alpha, beta, allele_count, total_mutations])
        ratio = beta / alpha
        log_likelihood = (
            ratio + 0.5 * (np.log(ratio) + np.log(2.0 / math.pi)) +
            allele_count * np.log(total_mutations / alpha) +
            (1 + 2 * allele_count) / 4 * (np.log(beta / (beta + 2 * total_mutations))) +
            compute_log_bessel(
                0.5 + allele_count, np.sqrt(beta * (beta + 2 * total_mutations)) / alpha
            )
        )
        return log_likelihood

    @staticmethod
    def compute_pdf(alpha, beta, allele_count, total_mutations, s, p_n=None):
        if(p_n is None):
            p_n = FitnessEstimator.log_likelihood_analytical(
                alpha, beta, allele_count, total_mutations
            )
        lam = total_mutations / s
        poisson_part = allele_count * math.log(lam) - lam
        inv_gaussian_part = 0.5 * math.log(beta / (2 * np.pi * math.pow(s, 3))) + (
            -beta * math.pow(s - alpha, 2) / (2 * alpha * alpha * s)
        )
        return math.exp(poisson_part + inv_gaussian_part - p_n)

    @staticmethod
    def compute_posterior(alpha, beta, allele_count, total_mutations, delta=0.0005, N_steps=5000):
        p_n = FitnessEstimator.log_likelihood_analytical(alpha, beta, allele_count, total_mutations)
        total_prob = 0
        posterior_pdf = np.zeros(N_steps)
        for i in range(N_steps):
            s = delta * (i + 1)
            posterior_pdf[i] = FitnessEstimator.compute_pdf(
                alpha, beta, allele_count, total_mutations, s, p_n
            )
            total_prob += delta * posterior_pdf[i]
        # print(str(total_prob) + ' ' + str(allele_count) + ' ' + str(total_mutations))
        '''
        plt.plot(delta * np.array(range(N)), np.log(posterior_pdf), 'o')
        plt.show()
        pdb.set_trace()
        '''
        return posterior_pdf

    @staticmethod
    def compute_mean(alpha, beta, allele_count, total_mutations, delta=0.0005):
        posterior_pdf = FitnessEstimator.compute_posterior(
            alpha, beta, allele_count, total_mutations, delta
        )
        total = 0
        for i in range(len(posterior_pdf)):
            s = delta * (i + 1)
            total += posterior_pdf[i] * delta * s
        return total

    # gene_df: a dataframe with columns of allele_count and mutation_rate
    def __init__(self, gene_df):
        # self._alpha, self._beta = self.estimate_prior(gene_df)
        self._alpha, self._beta = 0.05, 0.008
        # FitnessEstimator.plot_prior_and_posterior(0.05, 0.008, gene_df)

    @staticmethod
    def plot_prior_and_posterior(alpha, beta, gene_df):
        N_steps = 5000
        x_vals = np.linspace(0, 2.5, N_steps)
        prior_pdf = [scipy.stats.invgauss.pdf(x, alpha, 0, beta) for x in x_vals]
        posterior_pdf = np.zeros(N_steps)
        gene_count = 0
        for index, gene in gene_df.iterrows():
            posterior_pdf += FitnessEstimator.compute_posterior(
                alpha, beta, gene['allele_count'], gene['total_mutations'], N_steps=N_steps
            )
            gene_count += 1
            if(gene_count % 100 == 0):
                print(gene_count)
            if(gene_count > 1000):
                break
        posterior_pdf = posterior_pdf / gene_count
        plotting.plot_prior_posterior(x_vals, prior_pdf, posterior_pdf)

    def estimate_prior(self, gene_df):
        N_alpha = 40
        N_beta = 40
        alpha_vals = np.exp(np.linspace(np.log(.001), np.log(2), N_alpha))
        beta_vals = np.exp(np.linspace(np.log(.0001), np.log(2), N_beta))
        results = np.zeros((N_alpha, N_beta))
        for alpha_index in range(N_alpha):
            print(alpha_index)
            alpha = alpha_vals[alpha_index]
            for beta_index in range(N_beta):
                beta = beta_vals[beta_index]
                results[alpha_index, beta_index] = (
                    self.compute_log_likelihood(gene_df, alpha, beta)
                )
        np.save('Results/prior_results' + str(N_alpha) + '.npy', results)
        alpha_array, beta_array = np.meshgrid(alpha_vals, beta_vals, indexing='ij')
        plot_figure = True
        if(plot_figure):
            plotting.plot_prior_landscape(alpha_array, beta_array, results)
        max_index = np.ndarray.argmax(results)
        (ml_alpha_index, ml_beta_index) = np.unravel_index(max_index, results.shape)
        return alpha_vals[ml_alpha_index], beta_vals[ml_beta_index]

    def compute_log_likelihood(self, gene_df, alpha, beta):
        log_likelihood = 0
        number_genes = len(gene_df.index)
        for index, gene in gene_df.iterrows():
            log_likelihood += 1.0 / number_genes * FitnessEstimator.log_likelihood_analytical(
                alpha, beta, gene['allele_count'], gene['total_mutations']
            )
        return log_likelihood

    def calculate_posterior_mean(self, gene_df):
        gene_df['s_mean'] = -1
        gene_count = 0
        for index, gene in gene_df.iterrows():
            gene_df.loc[index, 's_mean'] = FitnessEstimator.compute_mean(
                self._alpha, self._beta, gene['allele_count'], gene['total_mutations']
            )
            gene_count += 1
            if(gene_count % 100 == 0):
                print(gene_count)
            if(gene_count > 1000):
                break
        pdb.set_trace()
