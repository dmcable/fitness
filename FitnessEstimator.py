from scipy.special import kv as bessel
import numpy as np
import pdb
import math
import scipy.stats
import plotting
import pandas as pd


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


def inv_gauss_cdf(alpha, beta, s_het):
    if(s_het < 1e-8):
        return 0
    if(beta / alpha > 300):
        return 1
    cdf = (
        scipy.stats.norm.cdf(math.sqrt(beta / s_het) * (s_het / alpha - 1)) +
        math.exp(2 * beta / alpha) *
        scipy.stats.norm.cdf(-math.sqrt(beta / s_het) * (s_het / alpha + 1))
    )
    return cdf


class FitnessEstimator(object):
    # gene_df: a dataframe with columns of allele_count and mutation_rate
    def __init__(self, gene_df):
        self._N_steps = 20000
        self._delta = 0.0001
        # self._alpha, self._beta = self.estimate_prior(gene_df, method='iterative', grid_size=40)
        self._alpha, self._beta = .027, .0027

    def estimate_prior(self, gene_df, method='empirical', grid_size=10):
        assert(method in ['empirical', 'iterative'])
        N_alpha = grid_size
        N_beta = grid_size
        alpha_vals = np.exp(np.linspace(np.log(.001), np.log(2), N_alpha))
        beta_vals = np.exp(np.linspace(np.log(.0001), np.log(2), N_beta))
        results = np.zeros((N_alpha, N_beta))
        for alpha_index in range(N_alpha):
            alpha = alpha_vals[alpha_index]
            for beta_index in range(N_beta):
                beta = beta_vals[beta_index]
                if(method == 'empirical'):
                    results[alpha_index, beta_index] = (
                        self.compute_log_likelihood(gene_df, alpha, beta)
                    )
                elif(method == 'iterative'):
                    self._alpha = alpha
                    self._beta = beta
                    results[alpha_index, beta_index] = self.prior_posterior_distance(gene_df)
        np.save('results/prior_results_simulated_' + str(N_alpha) + '_' + method + '.npy', results)
        alpha_array, beta_array = np.meshgrid(alpha_vals, beta_vals, indexing='ij')
        plot_figure = False
        if(plot_figure):
            plotting.plot_prior_landscape(alpha_array, beta_array, results)
        if(method == 'empirical'):
            best_index = np.ndarray.argmax(results)
        else:
            best_index = np.ndarray.argmin(results)
        (ml_alpha_index, ml_beta_index) = np.unravel_index(best_index, results.shape)
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
        for index, gene in gene_df.iterrows():
            try:
                gene_df.loc[index, 's_mean'] = FitnessEstimator.compute_mean(
                    self._alpha, self._beta, gene['allele_count'], gene['total_mutations']
                )
            except ValueError:
                gene_df.loc[index, 's_mean'] = 0

    def calculate_posterior_mean_strat(self, gene_df):
        strat_names = ['0', '1', '2', '3']
        for name in strat_names:
            gene_df['s_' + name] = -1
        for index, gene in gene_df.iterrows():
            for name in strat_names:
                try:
                    gene_df.loc[index, 's_' + name] = FitnessEstimator.compute_mean(
                        self._alpha, self._beta, gene['ac_' + name], gene['total_mutations']
                    )
                except ValueError:
                    gene_df.loc[index, 's_' + name] = 0

    def compare_priors(self, gene_df):
        prior_params = [
            (0.027, 0.0027), (0.05, 0.008), (0.011, 0.0038), (2.0, 0.0038), (0.016, 0.0038)
        ]
        labels = ['iterative best', 'paper', 'empirical min', 'empirical max', 'empirical best']
        s_est = pd.DataFrame()
        self.calculate_naive_mean(gene_df)
        s_est['naive'] = gene_df['s_naive']
        for i in range(len(labels)):
            self._alpha, self._beta = prior_params[i]
            self.calculate_posterior_mean(gene_df)
            s_est[labels[i]] = gene_df['s_mean']
        percent_errors = pd.DataFrame()
        for label in (labels + ['naive']):
            if label is not 'iterative best':
                percent_errors[label] = (
                    abs((s_est[label] - s_est['iterative best'])) / s_est['iterative best']
                )

    def calculate_naive_mean(self, gene_df):
        gene_df['s_naive'] = gene_df['total_mutations']/gene_df['allele_count']
        gene_df.replace(np.inf, 1, inplace=True)

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

    @staticmethod
    def load_results(grid_size, method):
        N_alpha = grid_size  # grid_size
        N_beta = grid_size  # grid_size
        alpha_vals = np.exp(np.linspace(np.log(.001), np.log(2), N_alpha))
        beta_vals = np.exp(np.linspace(np.log(.0001), np.log(2), N_beta))
        file_name = 'results/prior_results_simulated_' + str(N_alpha) + '_' + method + '.npy'
        results = np.load(file_name)
        if(method == 'empirical'):
            best_index = np.ndarray.argmax(results)
        else:
            best_index = np.ndarray.argmin(results)
        (ml_alpha_index, ml_beta_index) = np.unravel_index(best_index, results.shape)
        (best_alpha, best_beta) = alpha_vals[ml_alpha_index], beta_vals[ml_beta_index]
        return best_alpha, best_beta, results, alpha_vals, beta_vals

    def get_aggregate_posterior_cdf(self, gene_df):
        posterior_pdf = np.zeros(self._N_steps)
        for index, gene in gene_df.iterrows():
            posterior_pdf += FitnessEstimator.compute_posterior(
                self._alpha, self._beta, gene['allele_count'], gene['total_mutations'],
                N_steps=self._N_steps, delta=self._delta
            )
        posterior_pdf = posterior_pdf / len(gene_df)
        posterior_cdf = np.cumsum(posterior_pdf) * self._delta
        return posterior_cdf, posterior_pdf

    def fit_inverse_gaussian(self, pdf):
        first_moment = 0
        second_moment = 0
        for i in range(self._N_steps):
            s_val = self._delta * i
            first_moment += self._delta * pdf[i] * s_val
            second_moment += self._delta * pdf[i] * (s_val ** 2)
        norm_const = np.sum(pdf * self._delta)
        first_moment /= norm_const
        second_moment /= norm_const
        alpha_est = first_moment
        variance_est = second_moment - (first_moment ** 2)
        beta_est = (alpha_est ** 3) / variance_est
        return alpha_est, beta_est

    def plot_prior_and_posterior(self, gene_df):
        s_vals = np.array([i * self._delta for i in range(self._N_steps)])
        prior_cdf = [inv_gauss_cdf(self._alpha, self._beta, s) for s in s_vals]
        posterior_cdf, _ = self.get_aggregate_posterior_cdf(gene_df)
        self.calculate_naive_mean(gene_df)
        sorted_s_naive = np.sort(gene_df['s_naive'].values)
        plotting.plot_prior_posterior(
            s_vals, prior_cdf, posterior_cdf, take_log=False, s_naive=sorted_s_naive
        )

    def sample_prior(self, N_samples):
        s_vals = np.array([i * self._delta for i in range(self._N_steps)])
        prior_cdf = [inv_gauss_cdf(self._alpha, self._beta, s) for s in s_vals]
        return self._delta * np.searchsorted(prior_cdf, np.random.rand(N_samples))

    def iterative_prior(self, gene_df):
        _, posterior_pdf = self.get_aggregate_posterior_cdf(gene_df)
        self._alpha, self._beta = self.fit_inverse_gaussian(posterior_pdf)
        print('alpha: ' + str(self._alpha))
        print('beta: ' + str(self._beta))
        self.iterative_prior(gene_df)

    def prior_posterior_distance(self, gene_df):
        s_vals = np.array([i * self._delta for i in range(self._N_steps)])
        prior_cdf = [inv_gauss_cdf(self._alpha, self._beta, s) for s in s_vals]
        posterior_cdf, _ = self.get_aggregate_posterior_cdf(gene_df)
        # np.sum(self._delta * np.abs(prior_cdf - posterior_cdf))
        indices_outside_support = np.where(np.array(prior_cdf) > 0.999)[0]
        if(len(indices_outside_support) > 0):
            support_length = min(np.where(np.array(prior_cdf) > 0.999)[0][0] * self._delta, 1.0)
        else:
            support_length = 1.0
        return np.log(self._delta * np.sum(np.square(prior_cdf - posterior_cdf)) / support_length)
