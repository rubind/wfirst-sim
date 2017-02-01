from numpy import *
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cPickle as pickle
import sys
import corner
from scipy.stats import scoreatpercentile
import scipy.ndimage.filters as scifilters
import glob

import os
wfirst_path = os.environ["WFIRST"]

sys.path.append(wfirst_path + "/scripts/pixel-level/")
from pixel_level_ETC import save_img
sys.path.append(wfirst_path + "/scripts/cosmo/")
import astro_functions
from scipy.interpolate import interp1d
import multiprocessing as mp




def get_cred(vals, axis = 0):
    return scoreatpercentile(vals, 15.8655253931, axis = axis), median(vals, axis = axis), scoreatpercentile(vals, 84.13447460685, axis = axis)

def print_array(A):
    for row in A:
        for val in row:
            print '{:4}'.format(val),
        print
    

def distance_modulus_to_fractional(modulus):
    return modulus * log(10.)/5.

def convert_ax_mu_to_d(ax_mu):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax_mu.get_ylim()
    ax_d.set_ylim(distance_modulus_to_fractional(y1), distance_modulus_to_fractional(y2))
    ax_d.figure.canvas.draw()

def target_FoM_scale(coarse_z, coarse_mu_errs, target_FoM):
    err_scales = arange(0.1, 3., 0.1)
    FoMs_found = []

    for err_scale in err_scales:
        tmp_Cmat = diag(err_scale * coarse_mu_errs)**2.
        print tmp_Cmat
        FoM = astro_functions.get_FoM(tmp_Cmat, coarse_z, zp = 0.3, addcmb = 1, verbose = True, shift_constraint = 0.002)[0]
        FoMs_found.append(FoM)

    FoMs_found = array(FoMs_found)
    print "scales ", err_scales, FoMs_found
    inds = where(1 - isnan(FoMs_found))

    FoMs_found = FoMs_found[inds]
    err_scales = err_scales[inds]


    FoMs_found = FoMs_found[::-1]
    err_scales = err_scales[::-1]

    print "scales ", err_scales, FoMs_found

    interp_scale = interp1d(FoMs_found[inds], err_scales[inds], kind = 'linear')(target_FoM)
    print "interpolated to", interp_scale, "for", target_FoM

    return interp_scale * coarse_mu_errs



def make_coarse_plot(bin_edges, bin_center_floats, bin_centers, mu_errs, pltsuffix, title):
    plt.figure()

    for i in range(len(bin_edges) - 1):
        #if bin_center_floats[i] > 0.1:
        plt.plot(1 + array(bin_edges[i:i+2]), [mu_errs[i]/(5./log(10.))]*2, color = (0.8, 0.8, 0.8), zorder = 1)

    startat = 0 #int(bin_center_floats[0] < 0.1)
    plt.plot(1 + array(bin_center_floats[startat:]), mu_errs[startat:]/(5./log(10.)), color = 'orange', zorder = 2)
    plt.plot(1 + array(bin_center_floats[startat:]), mu_errs[startat:]/(5./log(10.)), 'o', color = 'orange', zorder = 3)

    print "Fractional Distance: ", bin_centers, mu_errs/(5./log(10.))

    plt.xscale('log')
    plt.xlim(1,4)
    plt.ylim(0,0.025)
    plt.title(title, size = 8)
    plt.savefig(outputdir + "dist_vs_z_coarse" + pltsuffix + ".pdf")
    plt.close()


def binned_points(the_cov, z_set, target_FoM, bin_edges = [0., 0.1, 0.2, 0.55, 0.85, 1.15, 1.7], pltsuffix = "", title = ""):
    the_weight = linalg.inv(the_cov)
    
    bin_jacob = zeros([len(z_set), len(bin_edges) - 1], dtype=float64)
    for i in range(len(bin_edges) - 1):
        for j in range(len(z_set)):
            if (z_set[j] > bin_edges[i])*(z_set[j] <= bin_edges[i+1]):
                bin_jacob[j, i] += 1
    
    print "bin_jacob:"
    print_array(bin_jacob)
    wmat = dot(transpose(bin_jacob), dot(the_weight, bin_jacob))
    try:
        cmat = linalg.inv(wmat)
    except:
        cmat = linalg.inv(wmat + diag(ones(len(wmat), dtype=float64))*0.000001)

    print_array(cmat)
    print sqrt(diag(cmat))
    
    bin_centers = []
    bin_center_floats = []
    for i in range(len(bin_edges) - 1):
        bin_centers.append("%.3f" % mean(bin_edges[i:i+2]))
        bin_center_floats.append(mean(bin_edges[i:i+2]))

    new_errs = target_FoM_scale(bin_center_floats, sqrt(diag(cmat)), target_FoM = target_FoM)
    
    make_coarse_plot(bin_edges=bin_edges, bin_center_floats=bin_center_floats, bin_centers=bin_centers, mu_errs = sqrt(diag(cmat)), pltsuffix = pltsuffix, title = title)
    make_coarse_plot(bin_edges=bin_edges, bin_center_floats=bin_center_floats, bin_centers=bin_centers, mu_errs = new_errs, pltsuffix = pltsuffix + "_rescaled", title = title)

    




def plot_delta_vs_value(true_vals, samples, pltname, diff = 1):
    plt.figure()
    c1, c2, c3 = get_cred(samples)

    plt.errorbar(true_vals, c2 - true_vals*diff, fmt = '.', capsize = 0, yerr = 0.5*(c3 - c1))
    plt.savefig(outputdir + pltname)
    plt.close()


def make_corner(keys, pltname, fit_params):
    plt.figure()
    samples = []
    labels = []
    
    for key in keys:
        if fit_params.has_key(key):
            if len(fit_params[key].shape) == 1:
                samples.append(fit_params[key])
                labels.append(key)
            else:
                for i in range(len(fit_params[key][0])):
                    samples.append(fit_params[key][:,i])
                    labels.append(key + "_" + str(i))
    corner.corner(transpose(array(samples)), labels = labels)
    plt.savefig(outputdir + pltname)
    plt.close()

def make_multi_corner(keys, pltname, fit_params, stan_data):
    pdf = PdfPages(outputdir + pltname)
    for i in range(stan_data["nred"]):
        labels = []
        samples = []
        for key in keys:
            if fit_params.has_key(key):
                samples.append(fit_params[key][:,i])
                labels.append(key)
        corner.corner(transpose(array(samples)), labels = labels)
        pdf.savefig(plt.gcf())
    pdf.close()
    plt.close()


def systematics_histograms(other_inputs, fit_params, stan_data):
    plt.figure(figsize = (6, 3*stan_data["nsys"]))
    for i in range(stan_data["nsys"]):
        plt.subplot(stan_data["nsys"], 1, i+1)
        plt.hist(fit_params["dsys"][:,i])
        c1, c2, c3 = get_cred(fit_params["dsys"][:,i])

        plt.title("%s %.3f + %.3f - %.3f" % (other_inputs["jacobian_names"][i], c2, c3 - c2, c2 - c1))

    plt.savefig(outputdir + "delta_sys_dist.pdf", bbox_inches = 'tight')
    plt.close()


def mag_histogram(fit_params):
    plt.figure()
    plt.hist(median(fit_params["mags"], axis = 0))
    plt.title("RMS %.3f" % std(median(fit_params["mags"], axis = 0), ddof = 1))
    plt.savefig(outputdir + "mags.pdf")
    plt.close()

def plot_d_mu_d_sys(other_inputs, stan_data, fit_params, offset = 0):
    color_dict = {'Count-Rate Nonlinearity': 'b', 'MW Extinction Normalization': 'r', 'MW Extinction Zeropoint': 'orange', 'MW Extinction $R_V$': 'g',
                  'IG Extinction': 'k', 'Fundamental Calibration': 'm', " ": 'k'}

    labeled = []
    for i in range(stan_data["nsys"]):
        dydx_vect = []
        for j in range(stan_data["nred"]):
            dydx = cov(fit_params["dsys"][:,i], fit_params["mus"][:,j])[0,1]/std(fit_params["dsys"][:,i])
            dydx_vect.append(dydx)
        dydx_vect = array(dydx_vect)
        ind = color_dict.keys().index(other_inputs["jacobian_names"][i])

        plt.figure(1)
        plt.plot(other_inputs["redshifts"], dydx_vect + ind*offset, color = color_dict[other_inputs["jacobian_names"][i]], label = other_inputs["jacobian_names"][i]*(labeled.count(other_inputs["jacobian_names"][i]) == 0))
        plt.figure(2)
        plt.plot(other_inputs["redshifts"], scifilters.gaussian_filter1d(dydx_vect, sigma = 2, mode = 'nearest') + ind*offset,
                 color = color_dict[other_inputs["jacobian_names"][i]], label = other_inputs["jacobian_names"][i]*(labeled.count(other_inputs["jacobian_names"][i]) == 0))
        labeled.append(other_inputs["jacobian_names"][i])

    for i in (1,2):
        plt.figure(i)
        plt.legend(loc = 'best', fontsize = 9, frameon = False)
        plt.xlabel("Redshift")
        plt.ylabel("Effect of Systematic (mag)")
        plt.savefig(outputdir + "dmu_dsys" + "_smoothed"*(i == 2) + "_offset"*(offset != 0) + ".pdf", bbox_inches = 'tight')
        plt.close()


def plot_d_mu_d_extinction(other_inputs, stan_data, fit_params):
    plt.figure()
    for key in ["RV_coeff", "log_R_EBV_coeff", "log_EBV_star_coeff"]:
        if fit_params.has_key(key):
            for i in range(stan_data["ncoeff"]):
                dydx_vect = []
                for j in range(stan_data["nred"]):
                    try:
                        dydx = cov(fit_params[key][:,i], fit_params["mus"][:,j])[0,1]/std(fit_params[key][:,i])
                    except:
                        dydx = -1
                        

                    dydx_vect.append(dydx)
                plt.plot(other_inputs["redshifts"], dydx_vect, label = "Coeff %i" % i)

            plt.legend(loc = 'best', fontsize = 9, frameon = False)
            plt.xlabel("Redshift")
            plt.ylabel("Effect of " + key.replace("_coeff", "") + " (mag)")
            plt.savefig(outputdir + "dmu_d" + key + ".pdf", bbox_inches = 'tight')
            plt.close()

def plot_RV_by_z(other_inputs, fit_params):
    plt.figure()

    c1, c2, c3 = get_cred(fit_params["RV_by_red"])
    plt.fill_between(other_inputs["redshifts"], c1, c3, color = 'cyan')
    plt.plot(other_inputs["redshifts"], c2, color = 'k')
    plt.xlabel("Redshift")
    plt.savefig(outputdir + "RV_credible_by_z.pdf", bbox_inches = 'tight')
    plt.close()

def plot_d_mu_d_RV_by_z(other_inputs, stan_data, fit_params):

    for i in range(stan_data["nred"]):

        the_cmat = cov(fit_params["RV_by_red"][:,i], fit_params["mus"][:,i])
        C_RVmu = the_cmat[0,1]
        C_RV = the_cmat[0,0]

        plt.figure(1)
        plt.plot(other_inputs["redshifts"][i], C_RVmu/C_RV, 'o', color = 'b')
        plt.figure(2)
        plt.plot(other_inputs["redshifts"][i], C_RVmu/sqrt(C_RV), 'o', color = 'b')


    plt.figure(1)
    plt.xlabel("Redshift")
    plt.ylabel("d mu/d RV")
    plt.savefig(outputdir + "dmu_dRV_z.pdf", bbox_inches = 'tight')
    plt.close()

    plt.figure(2)
    plt.xlabel("Redshift")
    plt.ylabel("d mu/d RV")
    plt.savefig(outputdir + "dmu_from_RV_z.pdf", bbox_inches = 'tight')
    plt.close()

def get_dmu_from_param(stan_data, fit_params, param):
    color_dict = {'Count-Rate Nonlinearity': 'b', 'MW Extinction Normalization': 'r',
                  'MW Extinction Zeropoint': 'orange', 'MW Extinction $R_V$': 'g',
                  'IG Extinction': 'k', 'Fundamental Calibration': 'm', " ": 'k'}
    

    this_param = reshape(fit_params[param], [fit_params[param].shape[0], prod(fit_params[param].shape[1:])])
    
    ndim = len(this_param[0])

    print param, ndim, this_param.shape

    dmus = zeros([ndim, stan_data["nred"]], dtype=float64)
    
    print "HERE1", this_param.shape, len(this_param.shape)
    for i in range(stan_data["nred"]):
        for j in range(ndim):
            if len(this_param.shape) == 2:
                #print len(this_param[:,j]), len(fit_params["mus"][:,i])
                the_cmat = cov(this_param[:,j], fit_params["mus"][:,i])

            #print "the_cmat", the_cmat
            assert the_cmat.shape == (2,2), "Bad shape! " + str(the_cmat.shape)
            C_xy = the_cmat[0,1]
            C_xx = the_cmat[0,0]
            dmus[j, i] = C_xy/sqrt(C_xx)
    print "HERE2"
    return dmus

        

def plot_d_mu_d_all_by_z(other_inputs, stan_data, fit_params):

    plt.figure(figsize = (8, 12))

    plot_params = [
        ["dsys"],
        ["RV_coeff", "log_R_EBV_coeff"]]
        #["proj_star_variable_coeff", "log_R_proj_coeff"]]

    n_coeff = stan_data["ncoeff"]

    plot_titles = ["Systematics", "Host-Galaxy Extinction"]# + ["Eigenvectors"]

    plot_names = [other_inputs["jacobian_names"],
                  ["$R_V$"]*n_coeff + ["$\\tau_{E(B-V)}$"]*n_coeff]

    plot_colors = {'Count-Rate Nonlinearity': 'b', 'MW Extinction Normalization': 'r',
                   'MW Extinction Zeropoint': 'orange', 'MW Extinction $R_V$': 'g',
                   'IG Extinction': 'k', 'Fundamental Calibration': 'm', " ": 'k',
                   '$R_V$': 'b', '$\\tau_{E(B-V)}$': 'r'}


    ylim = [10000, -10000]

    for i in range(len(plot_params)):
        plt.subplot(len(plot_params), 1, i+1)
        running_count = 0
        labeled_items = []

        for j in range(len(plot_params[i])):
            dmu = get_dmu_from_param(stan_data, fit_params, plot_params[i][j])
            print "dmu ", plot_params[i][j], dmu.shape

            for k in range(len(dmu)):

                the_name = plot_names[i][running_count]
                the_color = plot_colors[the_name]
                if labeled_items.count(the_name):
                    the_name = ""
                else:
                    labeled_items.append(the_name)

                plt.plot(other_inputs["redshifts"], dmu[k] - mean(dmu[k]), color = the_color, label = the_name)
                running_count += 1

        plt.xlabel("Redshift")
        plt.ylabel("Effect on Distance Modulus (mag)")
        plt.title(plot_titles[i])
        plt.legend(loc = 'best', fontsize = 10)
        ylim[0] = min(plt.ylim()[0], ylim[0])
        ylim[1] = max(plt.ylim()[1], ylim[1])
    
    for i in range(len(plot_params)):
        plt.subplot(len(plot_params), 1, i+1)
        plt.ylim(ylim)

    plt.savefig("all_dmu.eps", bbox_inches = 'tight')
    plt.close()
    


def plot_EBV_err(fit_params):
    plt.figure()
    c1, c2, c3 = get_cred(fit_params["true_EBV"], axis = 0)
    plt.plot(c2, 0.5*(c3 - c1), '.')
    plt.savefig(outputdir + "sigma_EBV_vs_EBV.pdf")
    plt.close()

def plot_diag_mu_err(other_inputs, stan_data, fit_params):
    plt.figure()
    c1, c2, c3 = get_cred(fit_params["mus"], axis = 0)

    for i in range(len(c1)):
        nsne = sum(array(stan_data["zinds"]) == i)
        plt.plot(other_inputs["redshifts"][i], 0.5*(c3[i] - c1[i])*sqrt(nsne), 'o')
        plt.text(other_inputs["redshifts"][i], 0.5*(c3[i] - c1[i])*sqrt(nsne), str(nsne), size = 8)
    plt.savefig(outputdir + "sigma_indiv_mu_vs_z.pdf")
    plt.close()


def make_plots_from_dir(item):
    print "Loading pickle file...", item

    if item.count("/") > 0:
        outputdir = item[:item.rfind("/")] + "/"
    else:
        outputdir = ""
    
    (fit_params, stan_data, other_inputs) = pickle.load(open(item, 'rb'))


    for key in fit_params:
        print "fit_params:", key
    for key in other_inputs:
        print "other_inputs:", key

    if not other_inputs.has_key("jacobian_names"):
        other_inputs["jacobian_names"] = " "*stan_data["nsys"]
        other_inputs["redshifts"] = arange(stan_data["nred"])*0.1 + 0.05
        print "Hack, but fixed in later versions!"


    cmat = cov(transpose(fit_params["mus"]))
    save_img(cmat, outputdir + "cmat.fits")

    z_set = concatenate(([0.05], arange(len(cmat) - 1, dtype=float64)*0.05 + 0.125))
    print "z_set", z_set
    FoM = astro_functions.get_FoM(cmat, z_set, zp = 0.3, addcmb = 1, verbose = True, shift_constraint = 0.002)[0]
    print "Final_FoM", FoM
    binned_points(cmat, z_set, bin_edges = [0., 0.1, 0.2, 0.55, 0.85, 1.15, 1.7], target_FoM = FoM)



    plot_RV_by_z(other_inputs = other_inputs, fit_params = fit_params)
    mag_histogram(fit_params = fit_params)
    plot_diag_mu_err(other_inputs = other_inputs, stan_data = stan_data, fit_params = fit_params)
    plot_EBV_err(fit_params = fit_params)

    plot_d_mu_d_all_by_z(other_inputs, stan_data, fit_params)

    systematics_histograms(other_inputs = other_inputs, fit_params = fit_params, stan_data = stan_data)
    plot_d_mu_d_sys(other_inputs = other_inputs, stan_data = stan_data, fit_params = fit_params)
    plot_d_mu_d_sys(other_inputs = other_inputs, stan_data = stan_data, offset = 0.02, fit_params = fit_params)
    plot_d_mu_d_extinction(other_inputs = other_inputs, fit_params = fit_params, stan_data = stan_data)
    plot_d_mu_d_RV_by_z(other_inputs = other_inputs, stan_data = stan_data, fit_params = fit_params)

    plot_delta_vs_value(other_inputs["true_EBVs"], fit_params["true_EBV"], "EBV_resid.pdf")
    plot_delta_vs_value(other_inputs["true_mags"], fit_params["mags"], "mag_resid.pdf")
    if fit_params.has_key("true_dRV"):
        plot_delta_vs_value(other_inputs["true_RVs"] - 3.1, fit_params["true_dRV"], "RV_resid.pdf", diff = 0)
    plot_delta_vs_value(median(fit_params["mus"], axis = 0), fit_params["mus"], "mus_vs_mus.pdf")


    print "Corner plots..."
    make_corner(["RV_coeff", "log_EBV_star_coeff", "log_R_EBV_coeff"], "EBV_params.pdf", fit_params = fit_params)
    make_multi_corner(["mus", "RV_by_red", "log_EBV_star_by_red", "log_R_EBV_by_red"], "RV_mus.pdf", fit_params = fit_params, stan_data = stan_data)
    return None


########################################## Starts Here ##########################################

print "Usage: python ../plot_results.py pickle_hybrid.txt"
redo_old_results = 0


files_to_do = sys.argv[1:]
for i in range(len(files_to_do))[::-1]:
    item = files_to_do[i]

    if item.count("/") > 0:
        outputdir = item[:item.rfind("/")] + "/"
    else:
        outputdir = ""

    if redo_old_results:
        pass
    else:
        if glob.glob(outputdir + "cmat.fits") != []:
            print "No need to redo ", files_to_do[i]
            del files_to_do[i]
   

procs_to_use = 20

pool = mp.Pool(processes = procs_to_use)
pool.map(make_plots_from_dir, files_to_do)

print "Done!"

