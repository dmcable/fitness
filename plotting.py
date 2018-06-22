import matplotlib.pyplot as plt
import numpy as np


def plot_prior_landscape(alpha_array, beta_array, results):
    plt.figure()
    cp = plt.pcolor(np.log(alpha_array), np.log(beta_array), results)
    plt.colorbar(cp)
    plt.title('log likelihood')
    plt.xlabel('log alpha')
    plt.ylabel('log beta')
    plt.show()


def plot_prior_posterior(x_vals, prior_pdf, posterior_pdf):
    plt.figure()
    plt.plot(x_vals, np.log(prior_pdf))
    plt.hold(True)
    plt.plot(x_vals, np.log(posterior_pdf), 'r')
    plt.show()
