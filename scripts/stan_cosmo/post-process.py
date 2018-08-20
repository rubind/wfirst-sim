import cPickle as pickle
from numpy import *
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import time
import sys
from DavidsNM import save_img, eig
import multiprocessing as mp
from functools import partial

import sys
import os
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/synth_dataset/")

from synth_functions import NMAD, get_weights


def make_sim_data():
    stan_data = {}
    fit_params = {}

    nsamp = 2000
    nsne = 100

    stan_data["norms_0"] = array([1., 2.511886])
    stan_data["projs_0"] = array([[-1.], [1.]])
    
    fit_params["true_x"] = zeros([nsamp, nsne, 1], dtype=float64)

    fit_params["true_dM"] = zeros(nsne, dtype=float64)
    fit_params["true_dM"][0::2] += 0
    fit_params["true_dM"][1::2] += -1

    fit_params["dM"] = zeros([nsamp, nsne], dtype=float64)
    fit_params["dM"] += random.normal(size = nsne)*0.1 + fit_params["true_dM"]


    fit_params["true_x"][:,::2,:] = -1
    fit_params["true_x"][:,1::2,:] = 1


    fit_params["xstar"] = zeros([nsamp, 1], dtype=float64)

    stan_data["nsne"] = nsne
    return fit_params, stan_data

def filter_xy(x, y):
    x = array(x)
    y = array(y)

    inds = where(1 - isnan(x) - isnan(y) - isinf(x) - isinf(y))
    return x[inds], y[inds]


def get_dMs(i, pfl, dP, projs_0, nsne, xstar_by_SN, true_x, dM, norms_0_in_mags, NMAD_projs):

    print "samp ", i, time.asctime(), "dP", dP, pfl

    dMs = []
    for j in range(nsne):
        weights = get_weights(true_x[i][j] + xstar_by_SN[i][j], projs_0 = projs_0, dP = dP, NMAD_projs = NMAD_projs, clip_kernel = 0)

        dMs.append(   sum(weights*(dM[i][j] - norms_0_in_mags))
                      /sum(weights)   )

        """
        plt.hist(fit_params["dM"][i], bins = 50)
        plt.savefig("tmp.pdf")
        plt.close()
        aaaaa

        if j == 10:
            plt.plot(fit_params["dM"][i][j] - norms_0_in_mags, weights, '.')
            plt.axvline(dMs[-1])
            plt.savefig("tmp.pdf")
            plt.close()
            fff
        """
    return dMs
        
clip_kernel = 0

pool = mp.Pool(processes = 4)


plt.figure(1, figsize = (12, 16))

tstart = time.time()
dPs = 10.**linspace(-0, 0.5, 2.)

runs_total = len(dPs)*len(sys.argv[1:])*2000
run_counter = 1


for pfl in sys.argv[1:]:
    print "Loading...", time.asctime()
    [stan_data, fit_params] = pickle.load(open(pfl, 'rb'))
    stan_data = pickle.load(open("input.pickle", 'rb'))
    print time.asctime()

    """
    
    fit_params, stan_data = make_sim_data()
    """

    NMAD_projs = NMAD(stan_data["projs_0"], axis = 0)

    norms_0_in_mags = -2.5*log10(stan_data["norms_0"])
    print "norms_0_in_mags", norms_0_in_mags

    assert all(1 - isnan(norms_0_in_mags))

    """
    for i in range(stan_data["nx0"]):
        print i
        for j in range(i+1, stan_data["nx0"]):
            delta = (stan_data["projs_0"][i] - stan_data["projs_0"][j])/NMAD_projs
            delta = sqrt(dot(delta, delta))
            plt.plot(delta, abs(norms_0_in_mags[i] - norms_0_in_mags[j]), '.', color = 'b', alpha = 0.1)
    plt.savefig("dM_vs_dP.pdf")
    plt.close()
    """


    """
    fit_params, stan_data = make_sim_data()
    """
    
    print fit_params["true_x"].shape
    xerrs = NMAD(fit_params["true_x"], axis = 0)
    xerrs = median(xerrs, 0)
    plt.figure(2)
    plt.plot(xerrs/NMAD_projs)

    plt.figure(1)
    #NMAD(vals, axis = 0)
    

    nsamp = 2000

    plty = []
    pltHR = []
    pltHRNMAD = []


    for dP in dPs:
        
        # matrix [nEV, nredcoeff] xstar => [nsamp, nEV, nredcoeff];
        # vector [nredcoeff] redcoeff [nsne] => [nsamp, nsne, nredcoeff];
        xstar_by_SN = []# [nsamp, nsne, nEV]
        for j in range(nsamp):
            xstar_by_SN.append(dot(stan_data["redcoeff"], transpose(fit_params["xstar"][j])))
        xstar_by_SN = array(xstar_by_SN)
        print "xstar_by_SN", xstar_by_SN.shape
        
        """
        all_dMs = pool.map(partial(get_dMs, pfl = pfl, dP = dP, projs_0 = stan_data["projs_0"], nsne = stan_data["nsne"], xstar_by_SN = xstar_by_SN,
                                   true_x = fit_params["true_x"], dM = fit_params["dM"], norms_0_in_mags = norms_0_in_mags, NMAD_projs = NMAD_projs), range(nsamp))
        """
        all_dMs = []
        for j in range(nsamp):
            all_dMs.append(get_dMs(j, pfl = pfl, dP = dP, projs_0 = stan_data["projs_0"], nsne = stan_data["nsne"], xstar_by_SN = xstar_by_SN,
                                   true_x = fit_params["true_x"], dM = fit_params["dM"], norms_0_in_mags = norms_0_in_mags, NMAD_projs = NMAD_projs)
                           )
        all_dMs = array(all_dMs)

        all_dMs = ma.masked_array(all_dMs, mask = isnan(all_dMs))

        save_img(array(all_dMs), "all_dMs.fits")
        cmat = ma.cov(transpose(all_dMs))
        
        """
        evals, evecs, evecs_norm = eig(cmat)
        save_img(cmat, "cmat_dP=%.2f.fits" % dP)
        save_img(evecs_norm, "evecs_dP=%.2f.fits" % dP)
        """
        
        print cmat.shape

        print all_dMs.shape

        """
        plt.hist(all_dMs[0], bins = 50, alpha = 0.5)
        plt.hist(all_dMs[10], bins = 50, alpha = 0.5)
        plt.savefig("tmp.pdf")
        plt.close()
        ffff
        """


        stan_data["redshifts"] = array(stan_data["redshifts"])
        zbins = sort(unique(stan_data["redshifts"]))


        h_resids = ma.median(all_dMs, axis = 0)
        assert len(h_resids) == stan_data["nsne"]

        rms_by_SN = zeros(stan_data["nsne"], dtype=float64)
        nmad_by_SN = zeros(stan_data["nsne"], dtype=float64)

        for i in range(len(zbins)):
            print zbins[i]
            
            inds = where(stan_data["redshifts"] == zbins[i])
            rms_by_SN[inds] = ma.std(h_resids[inds])
            nmad_by_SN[inds] = 1.4826*ma.median(abs(h_resids[inds] - ma.median(h_resids[inds])))

        plt.figure(2)
        plt.plot(rms_by_SN)
        plt.plot(nmad_by_SN)
        plt.savefig("rms_by_SN.pdf")
        plt.close()

        cmat_nmad = cmat + diag(nmad_by_SN**2.)
        cmat = cmat + diag(rms_by_SN**2.)
        jmat = ones([stan_data["nsne"], 1], dtype=float64)
        try:
            pwmat = dot(transpose(jmat), dot(linalg.inv(cmat), jmat))
            pcmat = float(1./sqrt(pwmat))
        except:
            print "Couldn't invert!"
            pcmat = 0.


        plty.append(pcmat)
        pltHR.append(std(h_resids, ddof = 1))
        pltHRNMAD.append(NMAD(h_resids))
        

        cmat = array(cmat, dtype=float64)
        cmat_nmad = array(cmat_nmad, dtype=float64)

        print "zbins ", list(zbins)
        jmat = zeros([stan_data["nsne"], len(zbins)], dtype=float64)
        for i in range(len(zbins)):
            print zbins[i]
            inds = where(stan_data["redshifts"] == zbins[i])
            jmat[inds, i] += 1

        save_img(jmat, "jmat.fits")
        try:
            pwmat = dot(transpose(jmat), dot(linalg.inv(cmat), jmat))
            pcmat = linalg.inv(pwmat)
            pwmat_nmad = dot(transpose(jmat), dot(linalg.inv(cmat_nmad), jmat))
            pcmat_nmad = linalg.inv(pwmat_nmad)
        except:
            pcmat = array([[1]])
            print "Couldn't invert!"

        withzbins = zeros([len(pcmat) + 1, len(pcmat)], dtype=float64)
        withzbins[1:,] = pcmat
        withzbins[0] = zbins

        save_img(withzbins, "pcmat_" + str(dP) + ".fits")

        withzbins_nmad = zeros([len(pcmat) + 1, len(pcmat)], dtype=float64)
        withzbins_nmad[1:,] = pcmat_nmad
        withzbins_nmad[0] = zbins

        save_img(withzbins_nmad, "pcmat_nmad_" + str(dP) + ".fits")


    plt.subplot(2,1,1)


    
    plt.plot(filter_xy(dPs, plty)[0], filter_xy(dPs, plty)[1], 'o', label = pfl)
    plt.xscale('log')
    plt.xlabel("dP")

    plt.subplot(2,1,2)
    plt.ylabel("Implied HD RMS")

    #plt.plot(dPs*1.02, pltHR, '^', label = pfl, alpha = 0.1)
    plt.plot(filter_xy(dPs, pltHRNMAD)[0], filter_xy(dPs, pltHRNMAD)[1], 'o', label = pfl + " NMAD")
    try:
        plt.text(dPs[-1], pltHRNMAD[-1], pfl + " NMAD", ha = 'right', va = 'center', fontsize = 8)
    except:
        pass
    plt.xscale('log')
    plt.xlabel("dP")

for i in range(1):
    plt.subplot(2,1,i+1)
    plt.legend(loc = 'best', fontsize = 7)
    plt.ylim(0, plt.ylim()[1])
plt.savefig("RMS_vs_dP" + "_clipkernel"*clip_kernel + ".pdf")
plt.close()

plt.figure(2)

plt.savefig("xerrs.pdf")
plt.close()
