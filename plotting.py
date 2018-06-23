import matplotlib.pyplot as plt
import numpy as np
import pdb


def plot_prior_landscape(alpha_array, beta_array, results):
    plt.figure()
    cp = plt.pcolor(np.log(alpha_array), np.log(beta_array), results)
    plt.colorbar(cp)
    plt.title('log likelihood')
    plt.xlabel('log alpha')
    plt.ylabel('log beta')
    plt.show()


def plot_prior_posterior(s_vals, prior_vals, posterior_vals, take_log=False, s_naive=None):
    plt.figure()
    graphs = []
    labels = []
    if(take_log):
        graphs += plt.plot(s_vals, np.log(prior_vals))
        labels.append('prior pdf')
    else:
        graphs += plt.plot(s_vals, prior_vals)
        labels.append('prior cdf')
    plt.hold(True)
    if(take_log):
        graphs += plt.plot(s_vals, np.log(posterior_vals), 'r')
        labels.append('posterior pdf')
    else:
        graphs += plt.plot(s_vals, posterior_vals, 'r')
        labels.append('posterior cdf')
    if(s_naive is not None):
        graphs += plt.plot(s_naive, np.arange(len(s_naive)) / float(len(s_naive)), 'g')
        labels.append('naive s_het cdf')
    plt.title('plot of prior and posterior')
    plt.xlabel('s')
    plt.ylabel('pdf or cdf')
    plt.legend(graphs, labels)
    plt.xlim(0, 1)
    plt.show()
