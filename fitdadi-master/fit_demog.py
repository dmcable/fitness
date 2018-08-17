import dadi
import numpy
import matplotlib.pyplot as plt

def make_plot(popt, ns, pts_l, data):
    fs = two_epoch(popt, ns, pts_l[2])
    true_dist = list(data)
    true_dist = numpy.array(true_dist)
    true_dist[numpy.isnan(true_dist)] = 0
    true_dist /= numpy.sum(true_dist)
    model_dist = numpy.array(fs)
    model_dist /= numpy.sum(model_dist)

    plt.plot(numpy.log(true_dist))
    plt.hold(True)
    plt.plot(numpy.log(model_dist))
    plt.show()


def two_epoch(params, ns, pts):
    nu, T = params
    gamma = 0
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs


def fit_demog():
    # Load the data
    data = dadi.Spectrum.from_file('sfs/sfs_' + str('syn') + '.txt')
    ns = data.sample_sizes

    # These are the grid point settings will use for extrapolation.
    pts_l = [ns, int(ns * 1.2), int(ns * 1.4)]

    # The Demographics1D and Demographics2D modules contain a few simple models,
    # mostly as examples. We could use one of those.
    # Instead, we'll work with our custom model
    func = two_epoch
    # Now let's optimize parameters for this model.

    # The upper_bound and lower_bound lists are for use in optimization.
    # Occasionally the optimizer will try wacky parameter values. We in particular
    # want to exclude values with very long times, very small population sizes, or
    # very high migration rates, as they will take a long time to evaluate.
    # Parameters are: (nu1F, nu2B, nu2F, m, Tp, T)
    upper_bound = [200, 5]
    lower_bound = [0.2, 0.005]

    # truth = [2, 0.05]

    p0 = [2, 0.05]
    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Perturb our parameters before optimization. This does so by taking each
    # parameter a up to a factor of two up or down.
    p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                                  lower_bound=lower_bound)
    # Do the optimization. By default we assume that theta is a free parameter,
    # since it's trivial to find given the other parameters. If you want to fix
    # theta, add a multinom=False to the call.
    # The maxiter argument restricts how long the optimizer will run. For real
    # runs, you will want to set this value higher (at least 10), to encourage
    # better convergence. You will also want to run optimization several times
    # using multiple sets of intial parameters, to be confident you've actually
    # found the true maximum likelihood parameters.
    print('Beginning optimization ************************************************')
    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                       lower_bound=lower_bound,
                                       upper_bound=upper_bound,
                                       verbose=len(p0), maxiter=3)
    numpy.save('popt_' + str(ns) + '.npy', popt)
    if(False):
        make_plot(popt, ns, pts_l, data)
    # Calculate the best-fit model AFS.
    model = func_ex(popt, ns, pts_l)
    # Likelihood of the data given the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, data)
    print('Maximum log composite likelihood: {0}'.format(ll_model))
    # The optimal value of theta given the model.
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    print('Optimal value of theta: {0}'.format(theta))

    numpy.save('theta_' + str(ns) + '.npy', theta)
    return popt, theta

if __name__ == '__main__':
    fit_demog()
