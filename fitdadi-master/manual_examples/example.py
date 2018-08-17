#! /usr/bin/env python

import numpy
import dadi
import Selection
import pdb
import time

#the demographic function we will be using for our analysis. This describes
#a two epoch model with one historical size change.
def two_epoch(params, ns, pts):
    nu, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T, nu)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

#set demographic parameters and theta. this is usually inferred from
#synonymous sites

demog_params = np.load('popt.npy').tolist() # [2, 0.05]
theta_ns = 4000.


#the demography + selection function. single size change and selection.
def two_epoch_sel(params, ns, pts): #
    nu, T, gamma = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs



#pdb.set_trace()
#integrate over a range of gammas without defining breaks

ns = numpy.array([249])
pts_l = [60, 80, 100]
spectra = Selection.spectra(demog_params, ns, two_epoch_sel, pts_l=pts_l,int_bounds=(1e-5,500), Npts=30, echo=True,mp=False)

#n_peoples = [100, 200, 400, 800, 1600, 3200, 6400]
#results = []

# for n_people in n_peoples:
#     ns = numpy.array([n_people])
#     base_pts = int(n_people * 1.5)
#     pts_l = [base_pts, int(base_pts * 1.2), int(base_pts * 1.4)]
#     start = time.time()
#     spectra = Selection.spectra(demog_params, ns, two_epoch_sel, pts_l=pts_l,int_bounds=(0.1, 5), Npts=5, echo=True, mp=True)
#     end = time.time()
#     time_elapsed = end - start
#     print('Runtime: ' + str(time_elapsed))
#     results.append(time_elapsed)
# pdb.set_trace()
#integrate over a range of gammas and set breaks
# pts_l = [600, 800, 1000]
# int_breaks = [1e-4, 0.1, 1, 100, 500]
# spectra = Selection.spectra(demog_params, ns, two_epoch_sel, pts_l=pts_l,int_breaks=int_breaks, Npts=300,echo=True, mp=True)


#the spectra can be pickled for usage later. This is especially convenient
#if the process of generating the spectra takes a long time.
#import pickle
#pickle.dump(spectra, open("example_spectra.sp","w"))

#load sample data
data = dadi.Spectrum.from_file('example.sfs')

#fit a DFE to the data
#initial guess and bounds
sel_params = [0.2, 1000.]
lower_bound = [1e-3, 1e-2]
upper_bound = [1, 50000.]
p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,
                              upper_bound=upper_bound)
pdb.set_trace()
popt = Selection.optimize_log(p0, data, spectra.integrate, Selection.gamma_dist,
                              theta_ns, lower_bound=lower_bound,
                              upper_bound=upper_bound, verbose=len(sel_params),
                              maxiter=30)

#get expected SFS for MLE
model_sfs = spectra.integrate(popt[1], Selection.gamma_dist, theta_ns)

#one possible characterization of the neutral+gamma DFE
def neugamma(mgamma, pneu, alpha, beta):
    mgamma=-mgamma
    #assume anything with gamma<1e-4 is neutral
    if (0 <= mgamma) and (mgamma < 1e-4):
        return pneu/(1e-4) + (1-pneu)*Selection.gamma_dist(-mgamma, alpha,
                                                           beta)
    else:
        return Selection.gamma_dist(-mgamma, alpha, beta) * (1-pneu)

#vectorize custom DFE
neugamma_vec = numpy.frompyfunc(neugamma, 4, 1)

sel_params = [0.2, 0.2, 1000.]
lower_bound = [1e-3, 1e-3, 1e-2]
upper_bound = [1, 1, 50000.]
p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,
                              upper_bound=upper_bound)
pdb.set_trace()
popt = Selection.optimize_log(p0, data, spectra.integrate, neugamma_vec,
                              theta_ns, lower_bound=lower_bound,
                              upper_bound=upper_bound, verbose=len(sel_params),
                              maxiter=30)

#another possible characterization of the neutral+gamma
def neugamma(mgamma, pneu, pgamma, alpha, beta):
    mgamma=-mgamma
    #assume anything with gamma<1e-4 is neutral
    if (0 <= mgamma) and (mgamma < 1e-4):
        return pneu/(1e-4) + pgamma*Selection.gamma_dist(-mgamma, alpha, beta)
    else:
        return Selection.gamma_dist(-mgamma, alpha, beta) * pgamma

#define a constraint function. At MLE, consfunc = 0
def consfunc(x, *args):
    return 1-sum(x[0:-2])

neugamma_vec = numpy.frompyfunc(neugamma, 5, 1)

sel_params = [0.2, 0.8, 0.2, 1000.]
lower_bound = [1e-3, 1e-3, 1e-3, 1e-2]
upper_bound = [1, 1, 1, 50000.]

p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,
                              upper_bound=upper_bound)
popt = Selection.optimize_cons(p0, data, spectra.integrate, neugamma_vec,
                               theta_ns, lower_bound=lower_bound,
                               upper_bound=upper_bound, verbose=len(sel_params),
                               maxiter=30, constraint=consfunc)
