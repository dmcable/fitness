from scipy.special import kv as bessel
import numpy as np
import matplotlib.pyplot as plt
import pdb
import math
import scipy.stats


def log_bessel_approx(n, x):
    return 1/2 * (np.log(np.pi) - np.log(2 * n)) - n * (np.log(np.e * x) - np.log(2 * n))


def compute_log_bessel(n, x):
    bessel_val = bessel(n, x)
    if(bessel_val == 0.0):
        return -np.inf
    log_bessel_val = np.log(bessel_val)
    if(np.isnan(log_bessel_val) or np.isinf(log_bessel_val)):
        raise ValueError
        log_bessel_val = log_bessel_approx(n, x)
    return log_bessel_val


def plot_bessel_approx():
    n = 1000
    x_vals = np.linspace(1, 10000, 1000)
    y_1 = [np.log(bessel(n, x_vals[i])) for i in range(1000)]
    y_2 = [log_bessel_approx(n, x_vals[i]) for i in range(1000)]
    plt.figure()
    plt.plot(x_vals, y_1, 'r')
    plt.hold(True)
    plt.plot(x_vals, y_2, 'b')
    plt.show()


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
            compute_log_bessel(0.5 + allele_count, np.sqrt(beta * (beta + 2 * total_mutations)) / alpha)
        )
        return log_likelihood

    @staticmethod
    def likelihood_analytical(alpha, beta, allele_count, total_mutations):
        assert_floats([alpha, beta, allele_count, total_mutations])
        ratio = beta / alpha
        probability = (
            np.exp(ratio) * math.sqrt(2 * ratio / math.pi) *
            math.pow(total_mutations / alpha, allele_count) *
            math.pow((beta / (beta + 2 * total_mutations)), (1 + 2 * allele_count) / 4) *
            bessel(0.5 + allele_count, math.sqrt(beta * (beta + 2 * total_mutations)) / alpha)
        )
        return probability

    @staticmethod
    def compute_posterior(alpha, beta, allele_count, total_mutations):
        p_n = FitnessEstimator.log_likelihood_analytical(alpha, beta, allele_count, total_mutations)
        delta = 0.0005
        N = 5000
        total_prob = 0
        posterior_pdf = np.zeros(N)
        for i in range(N):
            s = delta * (i + 1)
            lam = total_mutations / s
            poisson_part = allele_count * math.log(lam) - lam
            inv_gaussian_part = 0.5 * math.log(beta / (2 * np.pi * math.pow(s, 3))) + (
                -beta * math.pow(s - alpha, 2) / (2 * alpha * alpha * s)
            )
            posterior_pdf[i] = math.exp(poisson_part + inv_gaussian_part - p_n)
            total_prob += delta * posterior_pdf[i]
        #print(str(total_prob) + ' ' + str(allele_count) + ' ' + str(total_mutations))
        '''
        plt.plot(delta * np.array(range(N)), np.log(posterior_pdf), 'o')
        plt.show()
        pdb.set_trace()
        '''
        return posterior_pdf

    @staticmethod
    def log_likelihood_integral(alpha, beta, allele_count, total_mutations):
        assert_floats([alpha, beta, allele_count, total_mutations])
        delta = 0.0001
        N = 10000
        total_prob = 0
        results = [0 for i in range(N)]
        for i in range(N):
            s = delta * (i + 1)
            lam = total_mutations / s
            poisson_part = math.pow(lam, allele_count) * np.exp(-lam)
            inv_gaussian_part = math.sqrt(beta / (2 * np.pi * math.pow(s, 3))) * np.exp(
                -beta * math.pow(s - alpha, 2) / (2 * alpha * alpha * s)
            )
            results[i] = poisson_part * inv_gaussian_part
            total_prob += delta * poisson_part * inv_gaussian_part
        '''
        plt.figure()
        plt.plot(results)
        plt.show()
        pdb.set_trace()
        '''
        return np.log(total_prob)


    # gene_df: a dataframe with columns of allele_count and mutation_rate
    def __init__(self, N_genomes, gene_df):
        self.N_genomes = N_genomes
        #self._alpha, self._beta = self.estimate_prior(gene_df)
        FitnessEstimator.plot_prior_and_posterior(0.05, 0.008, gene_df)

    @staticmethod
    def plot_prior_and_posterior(alpha, beta, gene_df):
        x_vals = np.linspace(0, 2.5, 5000)
        pdf = [scipy.stats.invgauss.pdf(x, alpha, 0, beta) for x in x_vals]
        FitnessEstimator.compute_posterior(alpha, beta, 20.0, 0.05)
        my_pdf = np.zeros(5000)
        gene_count = 0
        for index, gene in gene_df.iterrows():
            allele_count = gene['allele_count']
            if(allele_count < 60):
                total_mutations = 60700 * gene['mut_rate']
                my_pdf += (
                    FitnessEstimator.compute_posterior(alpha, beta, allele_count, total_mutations)
                )
            gene_count += 1
            if(gene_count % 100 == 0):
                print(gene_count)
            if(gene_count > 1000):
                break
        posterior = my_pdf / gene_count
        pdb.set_trace()
        plt.figure()
        plt.plot(x_vals, np.log(pdf))
        plt.hold(True)
        plt.plot(x_vals, np.log(posterior), 'r')
        plt.show()

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
        #np.save('results20.npy', results)
        alpha_array, beta_array = np.meshgrid(alpha_vals, beta_vals, indexing ='ij')
        plt.figure()
        cp = plt.pcolor(np.log(alpha_array), np.log(beta_array), results)
        plt.colorbar(cp)
        plt.title('log likelihood')
        plt.xlabel('log alpha')
        plt.ylabel('log beta')
        plt.show()
        pdb.set_trace()
        (ml_alpha_index, ml_beta_index) = np.unravel_index(np.ndarray.argmax(results), results.shape)
        (ml_alpha, ml_beta) = (alpha_vals[ml_alpha_index], beta_vals[ml_beta_index])
        return ml_alpha, ml_beta

    def compute_log_likelihood(self, gene_df, alpha, beta):
        log_likelihood = 0
        number_genes = len(gene_df.index)
        for index, gene in gene_df.iterrows():
            allele_count = gene['allele_count']
            if(allele_count < 60):
                total_mutations = self.N_genomes * gene['mut_rate']
                log_likelihood += 1.0 / number_genes * (
                    FitnessEstimator.log_likelihood_analytical(alpha, beta, allele_count, total_mutations)
                )
        return log_likelihood
