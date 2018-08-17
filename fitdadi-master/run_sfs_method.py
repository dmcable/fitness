import numpy
import dadi
import Selection
import pdb
from fit_demog import fit_demog
import json
import pickle


# the demography + selection function. single size change and selection.
def two_epoch_sel(params, ns, pts):
    nu, T, gamma = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs


def create_sfs_files(next_N, prefix, N_strata=4):
    orig_N = 402832
    for i in range(N_strata):
        my_strat = str(i)
        if(i == 4):
            my_strat = 'syn'
        sfs_orig = numpy.load(prefix + my_strat + '.npy')
        sfs = numpy.zeros(next_N)
        for i in range(len(sfs_orig)):
            if(sfs_orig[i] > 0):
                target = max(min(1.0 * i / orig_N, 1), 0)
                for j in range(int(sfs_orig[i])):
                    sfs[min(int(numpy.random.binomial(next_N, target)), next_N-1)] += 1
        f = open('sfs/sfs_' + my_strat + '.txt', "w")
        my_len = next_N
        f.write(str(my_len) + ' unfolded\n')
        for j in range(my_len):
            f.write(str(sfs[j]) + ' ')
        f.close()


def load_specra_and_demog(next_N):
    demog_loaded = True
    if(not demog_loaded):
        demog_params, theta_ns = fit_demog()
    else:
        demog_params = numpy.load('popt_' + str(next_N) + '.npy')
        theta_ns = numpy.load('theta_' + str(next_N) + '.npy')

    n_people = next_N
    ns = numpy.array([n_people])
    base_pts = int(n_people)
    pts_l = [base_pts, int(base_pts * 1.2), int(base_pts * 1.4)]
    gamma_bounds = (1e-5, 500)
    spectra_loaded = True
    if(not spectra_loaded):
        spectra = Selection.spectra(demog_params, ns, two_epoch_sel, pts_l=pts_l, int_bounds=gamma_bounds, Npts=30, echo=True, mp=True)
        pickle.dump(spectra, open('spectra.pkl', 'wb'))
    else:
        spectra = pickle.load(open('spectra.pkl', "rb"))
    return theta_ns, spectra


def fit_sel_params(next_N, prefix, theta_ns, integrate_fun, loaded_sfs=False, save=False):
    if(not loaded_sfs):
        create_sfs_files(next_N, prefix, 4)
    sel_params_opt = {}
    sel_params = [0.2, 1000.]
    lower_bound = [1e-3, 1e-2]
    upper_bound = [1, 50000.]
    p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,upper_bound=upper_bound)
    for i in range(4):
        my_strat = str(i)
        data = dadi.Spectrum.from_file('sfs/sfs_' + my_strat + '.txt')
        popt = Selection.optimize_log(p0, data, integrate_fun, Selection.gamma_dist,
                                      theta_ns, lower_bound=lower_bound,
                                      upper_bound=upper_bound, verbose=len(sel_params),
                                      maxiter=30)
        sel_params_opt[my_strat] = popt[1].tolist()
    if(save):
        with open('sel_params_' + str(next_N) + '.json', 'w') as fp:
            json.dump(sel_params_opt, fp)
    return sel_params_opt


def fit_gene_sel(next_N, gene_list, theta_ns, integrate_fun):
    N_quartiles = 4
    N_params = 2
    gene_params = numpy.zeros((len(gene_list), N_quartiles, N_params))
    for i in range(len(gene_list)):
        gene = gene_list[i]
        prefix = 'gene_sfs/sfs_' + gene + '_strat_'
        sel_params_opt = fit_sel_params(next_N, prefix, theta_ns, integrate_fun)
        for quart in range(N_quartiles):
            gene_params[i, quart, :] = sel_params_opt[str(quart)]
    numpy.save('gene_params', gene_params)
    print(gene_params)
    pdb.set_trace()


if __name__ == '__main__':
    created_sfs = False
    next_N = 10000
    if(not created_sfs):
        create_sfs_files(next_N, 'sfs/sfs_fullest_strat_', 5)
    theta_ns, spectra = load_specra_and_demog(next_N)
    sel_params = fit_sel_params(next_N, 'sfs/sfs_fullest_strat_', theta_ns, spectra.integrate, True, True)
    print(sel_params)
    pdb.set_trace()
    gene_list = ['TTN'] # ['BRPF3', 'MOCOS']
    fit_gene_sel(next_N, gene_list, theta_ns, spectra.integrate)
