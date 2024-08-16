import numpy as np
import pickle
import matplotlib.pyplot as plt
import extinction
from DavidsNM import save_img, miniLM_new
from scipy.interpolate import RectBivariateSpline, interp1d
import tqdm
import argparse


def load_data():
    SN_data = pickle.load(open(opts.pickle, 'rb'))

    for key in SN_data:
        print("SN_data", key)

    print("SN_observations[0]", SN_data["SN_observations"][0])

    print("SN_table", SN_data["SN_table"])

    made_plot = 0

    plt.figure(figsize = (12, 4))
    for i in range(len(SN_data["SN_observations"])):
        if len(SN_data["SN_observations"][i]["dates"]) > 40 and made_plot == 0:
            SNRs = SN_data["SN_observations"][i]["fluxes"]/SN_data["SN_observations"][i]["dfluxes"]

            if np.median(SNRs) < 8:
                made_plot = 1
                for filt in np.unique(SN_data["SN_observations"][i]["filts"]):
                    inds = np.where(SN_data["SN_observations"][i]["filts"] == filt)
                    plt.subplot(1,2,1)
                    plt.plot(SN_data["SN_observations"][i]["dates"][inds], SN_data["SN_observations"][i]["fluxes"][inds])
                    plt.subplot(1,2,2)
                    plt.plot(SN_data["SN_observations"][i]["dates"][inds], SN_data["SN_observations"][i]["true_fluxes"][inds])
    plt.savefig("SN_LC.pdf", bbox_inches = 'tight')
    plt.close()

    return SN_data

def rest_model_1D_to_2D(rest_1D, other_data):
    rest_2D = np.zeros([len(other_data["model_phase"]),
                        len(other_data["model_rest"])], dtype=np.float64)

    ind_1D = 0
    for i in range(len(other_data["model_phase"])):
        for j in range(len(other_data["model_rest"])):
            if (i == other_data["phase0_ind"]) and ((j == other_data["restB_ind"]) or (j == other_data["restV_ind"])):
                pass
            else:
                rest_2D[i, j] = rest_1D[ind_1D]
                ind_1D += 1
    assert ind_1D == len(rest_1D)
    return rest_2D

def bin_by_wave(x, y, sigy, bin_edges):
    bin_vals = []
    bin_sigvals = []

    weighty = 1./sigy**2.
    
    for i in range(len(bin_edges) - 1):
        inds = np.where((x >= bin_edges[i])*(x < bin_edges[i+1]))
        bin_vals.append(sum(weighty[inds]*y[inds])/sum(weighty[inds]))
        bin_sigvals.append(1./np.sqrt(sum(weighty[inds])))

    bin_vals = np.array(bin_vals)
    bin_sigvals = np.array(bin_sigvals)

    assert np.isclose(sum(1./bin_sigvals**2.), sum(weighty))
    return bin_vals, bin_sigvals
    

def get_stan_data(SN_data):
    other_data = dict(filt_names = ["R062", "Z087", "Y106", "J129", "H158", "F184", "K213", "g", "r", "i", "z"],
                      lambs_dict = dict(R062 = 6200., Z087 = 8700., Y106 = 10600., J129 = 12900., H158 = 15800., F184 = 18400., K213 = 21300., g = 4800., r = 6200., i = 7700., z = 8700.),
                      milliz_list = [], model_phase = np.linspace(-15., 45., 10), model_rest = np.exp(np.linspace(np.log(3000.), np.log(16000.), opts.model_res)))

    other_data["phase0_ind"] = np.argmin(np.abs(other_data["model_phase"]))
    other_data["restB_ind"] = np.argmin(np.abs(other_data["model_rest"] - 4400.))
    other_data["restV_ind"] = np.argmin(np.abs(other_data["model_rest"] - 5500.))

    other_data["wave_bin_edges"] = np.exp(np.linspace(np.log(SN_data["IFC_waves"][0]*0.9999),
                                                      np.log(SN_data["IFC_waves"][-1]*1.0001),
                                                      prism_bins+1))
                                          

    other_data["AB25"] = (0.10884806248e-10/SN_data["IFC_waves"]**2.)
    print("wave_bin_edges", other_data["wave_bin_edges"])
    other_data["wave_bin"] = bin_by_wave(x = SN_data["IFC_waves"], y = SN_data["IFC_waves"], sigy = np.ones(len(SN_data["IFC_waves"])), bin_edges = other_data["wave_bin_edges"])[0]
    print("wave_bin", other_data["wave_bin"])


    prism_calib_fns = {}
    prism_calib_nodes = ["R062", "Z087", "Y106", "J129", "H158", "F184"]
    x_vals = [other_data["lambs_dict"][filt] for filt in prism_calib_nodes]

    for filt in prism_calib_nodes:
        y_vals = np.zeros(len(x_vals), dtype=np.float64)
        y_vals[prism_calib_nodes.index(filt)] = 1.
        prism_calib_fns[filt] = interp1d(x_vals, y_vals, kind = 'quadratic')
        
    
    stan_data = dict(NSNe = 0,
                     NFilt = len(other_data["filt_names"]),
                     redshifts = [],
                     daymaxes = [],
                     rest_lambs = [],
                     true_phases = [],
                     color_law = [],
                     gray_disp = [],
                     mags = [],
                     dmags = [],
                     z_inds = [],
                     d_mag_d_filt = [],
                     filt_inds = [])

    assert len(SN_data["SN_observations"]) == len(SN_data["SN_table"]["redshifts"])

    for SN_ind in range(len(SN_data["SN_observations"])):
        SN_LC = SN_data["SN_observations"][SN_ind]
        this_z = SN_data["SN_table"]["redshifts"][SN_ind]
        this_rest_lambs = np.array([other_data["lambs_dict"][item]/(1. + this_z) for item in SN_LC["filts"]])

        good_inds = np.where((this_rest_lambs > other_data["model_rest"][0])*
                             (this_rest_lambs < other_data["model_rest"][-1])*
                             (1 - np.isnan(SN_LC["true_fluxes"]))*
                             (SN_LC["true_fluxes"] > 0)
                             )

        this_rest_lambs = this_rest_lambs[good_inds]
        
        SNRs = SN_LC["true_fluxes"][good_inds] / (SN_LC["dfluxes"][good_inds])
        NPts = len(SN_LC["dates"][good_inds])


        if opts.SNRMax == 0:
            good_LC = (NPts > 10)*(np.sqrt(sum(SNRs**2.)) > 40)
        else:
            unique_filts = np.unique(SN_LC["filts"][good_inds])
            SNR_maxes = []
            for unique_filt in unique_filts:
                SNR_rest_lambs = np.array([other_data["lambs_dict"][item]/(1. + this_z) for item in SN_LC["filts"]])
                SNR_inds = np.where((SNR_rest_lambs > other_data["model_rest"][0])*
                                    (SNR_rest_lambs < other_data["model_rest"][-1])*
                                    (1 - np.isnan(SN_LC["true_fluxes"]))*
                                    (SN_LC["true_fluxes"] > 0)*(SN_LC["filts"] == unique_filt)
                                    )
                SNR_tmps = SN_LC["true_fluxes"][SNR_inds] / (SN_LC["dfluxes"][SNR_inds])
                SNR_maxes.append(np.max(SNR_tmps))
            SNR_maxes.sort()
            SNR_maxes = SNR_maxes[::-1]

            if len(SNR_maxes) > 2:
                good_LC = (SNR_maxes[0] >= 10)*(SNR_maxes[2] >= 5)
            else:
                good_LC = 0
                
        if good_LC:
            stan_data["mags"].append(-2.5*np.log10(SN_LC["true_fluxes"][good_inds]))
            stan_data["dmags"].append((2.5/np.log(10.)) * np.abs(SN_LC["dfluxes"][good_inds]/SN_LC["true_fluxes"][good_inds]))
            stan_data["d_mag_d_filt"].append(np.zeros([NPts, stan_data["NFilt"]], dtype=np.float64))

            this_daymax = SN_data["SN_table"]["daymaxes"][SN_ind]
            
            stan_data["redshifts"].append(this_z)
            stan_data["daymaxes"].append(this_daymax)
            

            for filt in np.unique(SN_LC["filts"][good_inds]):
                assert other_data["filt_names"].count(filt) == 1
                filt_mask = SN_LC["filts"][good_inds] == filt
                stan_data["d_mag_d_filt"][-1][:, other_data["filt_names"].index(filt)] = filt_mask

                
            milliz = int(np.around(1000.*this_z))
            if other_data["milliz_list"].count(milliz) == 0:
                other_data["milliz_list"].append(milliz)
            
            stan_data["z_inds"].append(other_data["milliz_list"].index(milliz))        
            stan_data["filt_inds"].append([other_data["filt_names"].index(item) for item in SN_LC["filts"][good_inds]])

            

            
            stan_data["rest_lambs"].append(this_rest_lambs)
            stan_data["true_phases"].append((SN_LC["dates"][good_inds] - this_daymax)/(1. + this_z))
            stan_data["color_law"].append(extinction.ccm89(this_rest_lambs, a_v = 3.1, r_v = 3.1))


            total_SNR_squared = 0.

            if len(SN_LC["IFS_dates"]) > 1:
                
                for IFS_ind in range(len(SN_LC["IFS_dates"])):
                    #print(SN_LC)
                    
                    this_rest_lambs = other_data["wave_bin"]/(1. + this_z)
                    bin_flux, bin_sig_flux = bin_by_wave(x = SN_data["IFC_waves"],
                                                         y = SN_LC["IFS_true_fluxes"][IFS_ind]/other_data["AB25"],
                                                         sigy = SN_LC["IFS_dfluxes"][IFS_ind]/other_data["AB25"],
                                                         bin_edges = other_data["wave_bin_edges"])

                    total_SNR_squared += sum((bin_flux/bin_sig_flux)**2.)

                    
                    good_inds = np.where((bin_flux > 0)
                                         *(this_rest_lambs > other_data["model_rest"][0])
                                         *(this_rest_lambs < other_data["model_rest"][-1]))

                    N_prism_waves = len(this_rest_lambs[good_inds])
                    frac_flux_err = bin_sig_flux[good_inds]/bin_flux[good_inds]

                    this_d_mag_d_filt = np.zeros([N_prism_waves, stan_data["NFilt"]], dtype=np.float64)

                    for filt in prism_calib_fns:
                        this_d_mag_d_filt[:, other_data["filt_names"].index(filt)] = prism_calib_fns[filt](other_data["wave_bin"][good_inds])
                                          
                    stan_data["d_mag_d_filt"][-1] = np.concatenate((stan_data["d_mag_d_filt"][-1], this_d_mag_d_filt))

                    assert N_prism_waves == len(frac_flux_err)
                    
                    stan_data["mags"][-1] = np.concatenate((
                        stan_data["mags"][-1],
                        -2.5*np.log10(bin_flux[good_inds])
                        ))
                    
                    stan_data["dmags"][-1] = np.concatenate((
                        stan_data["dmags"][-1],
                        (2.5/np.log(10.)) * np.abs(frac_flux_err)
                    ))

                    stan_data["rest_lambs"][-1] = np.concatenate((
                        stan_data["rest_lambs"][-1],
                        this_rest_lambs[good_inds]
                    ))
                    
                    stan_data["color_law"][-1] = np.concatenate((
                        stan_data["color_law"][-1],
                        extinction.ccm89(this_rest_lambs[good_inds], a_v = 3.1, r_v = 3.1)
                    ))
                    
                    stan_data["filt_inds"][-1] = np.concatenate((
                        stan_data["filt_inds"][-1], [-1]*N_prism_waves
                    ))

                    
                    stan_data["true_phases"][-1] = np.concatenate((
                        stan_data["true_phases"][-1], [(SN_LC["IFS_dates"][IFS_ind] - this_daymax)/(1. + this_z)]*N_prism_waves
                        ))
                    
                
            if opts.twins == 0:
                stan_data["gray_disp"].append(opts.gray_disp)
            else:
                this_gray_disp = np.sqrt(0.004138995017871641 + 12.627569849602615/(687.738490452654 + total_SNR_squared))
                stan_data["gray_disp"].append(this_gray_disp)

            
            stan_data["NSNe"] += 1



    
    stan_data["Nz"] = len(other_data["milliz_list"])
    for key in stan_data:
        try:
            stan_data[key][0][0]
            for i in range(len(stan_data[key])):
                stan_data[key][i] = np.array(stan_data[key][i])
        except:
            pass

    plt.figure(figsize = (12, 4))
    plt.subplot(1,2,1)
    plt.hist(np.concatenate(stan_data["rest_lambs"]), bins = 50)
    
    plt.subplot(1,2,2)
    plt.hist(np.concatenate(stan_data["true_phases"]), bins = 50)
    plt.savefig("lambs_and_phases.pdf", bbox_inches = 'tight')
    plt.close()

    stan_data["NCoeff"] = len(other_data["model_rest"])*len(other_data["model_phase"]) - 2
    stan_data["dmagdcoeff"] = []

    print("NCoeff", stan_data["NCoeff"])
    print("Computing matrix...")

    for SN_ind in tqdm.trange(stan_data["NSNe"]):
        Npts_this_SN = len(stan_data["mags"][SN_ind])
        
        this_dmagdcoeff = np.zeros([Npts_this_SN, stan_data["NCoeff"]], dtype=np.float64)

        for i in range(stan_data["NCoeff"]):
            rest_1D = np.zeros(stan_data["NCoeff"], dtype=np.float64)
            rest_1D[i] = 1.

            rest_2D = rest_model_1D_to_2D(rest_1D = rest_1D, other_data = other_data)
            ifn = RectBivariateSpline(other_data["model_phase"], other_data["model_rest"]/10000., rest_2D, kx = 2, ky = 2)

            for j in range(Npts_this_SN):
                this_dmagdcoeff[j, i] = ifn(stan_data["true_phases"][SN_ind][j], stan_data["rest_lambs"][SN_ind][j]/10000.)
            
            
        stan_data["dmagdcoeff"].append(this_dmagdcoeff)

    
    #print("stan_data", stan_data)
    print("other_data", other_data)

    plt.hist(stan_data["redshifts"], bins = np.arange(0, 2.51, 0.05))
    plt.savefig("redshifts_selected.pdf", bbox_inches = 'tight')
    plt.close()

    plt.plot(stan_data["redshifts"], stan_data["gray_disp"], '.')
    plt.savefig("gray_disp_vs_z.pdf", bbox_inches = 'tight')
    plt.close()
    

    f = open("redshifts_selected.txt", 'w')
    for redshift in np.sort(np.unique(stan_data["redshifts"])):
        f.write("%.3f  %i\n" % (redshift, sum(stan_data["redshifts"] == redshift)))
    f.close()
    
    
    return stan_data, other_data


def parseP(P, stan_data):
    parsed = {}
    ind = 0

    parsed["mu_bins"] = P[ind: ind+stan_data["Nz"]]
    ind += stan_data["Nz"]

    parsed["dZPs"] = P[ind: ind+stan_data["NFilt"]]
    ind += stan_data["NFilt"]

    parsed["coeff"] = P[ind: ind+stan_data["NCoeff"]]
    ind += stan_data["NCoeff"]

    return parsed

def unparseP(parsed):
    return np.concatenate((parsed["mu_bins"], parsed["dZPs"], parsed["coeff"]))

def residfn(P, wrapped_data):
    stan_data = wrapped_data[0]
    parsed = parseP(P, stan_data)
    
    resid_norm = []
    
    for i in tqdm.trange(stan_data["NSNe"]):
        if direct_inverse:
            cmat = np.diag(stan_data["dmags"][i]**2.)
            cmat += stan_data["gray_disp"][i]**2.
            cmat += np.outer(stan_data["color_law"][i]*0.5, stan_data["color_law"][i]*0.5) # +- 0.5 mag E(B-V)
        
            for filt in np.unique(stan_data["filt_inds"][i]):
                if filt != -1:
                    filt_mask = (stan_data["filt_inds"][i] == filt)
                    opt_mask = (stan_data["rest_lambs"][i] < 9000.)
                    cmat += np.outer(filt_mask*opt_mask*opts.color_scatter_opt, filt_mask*opt_mask*opts.color_scatter_opt)
                    cmat += np.outer(filt_mask*(1 - opt_mask)*opts.color_scatter_nir, filt_mask*(1 - opt_mask)*opts.color_scatter_nir)

            wmat = np.linalg.inv(cmat)
        else:
            # use Woodbury matrix identity

            Ainv = np.diag(stan_data["dmags"][i]**(-2.))
            unique_filt_inds = np.unique(stan_data["filt_inds"][i])
            Umat = np.zeros([len(stan_data["dmags"][i]), len(unique_filt_inds) + 2], dtype=np.float64)
            Umat[:, 0] = stan_data["gray_disp"][i]
            Umat[:, 1] = stan_data["color_law"][i]*0.5

            next_ind = 2
            for filt in unique_filt_inds:
                if filt != -1:
                    filt_mask = (stan_data["filt_inds"][i] == filt)
                    opt_mask = (stan_data["rest_lambs"][i] < 9000.)
                    Umat[:, next_ind] += filt_mask*opt_mask*opts.color_scatter_opt
                    Umat[:, next_ind] += filt_mask*(1 - opt_mask)*opts.color_scatter_nir
                    next_ind += 1
            Vmat = Umat.T

            identity_matrix = np.identity(len(Vmat), dtype=np.float64)
            middle_part = np.linalg.inv(identity_matrix + np.dot(Vmat, np.dot(Ainv, Umat)))
            wmat = Ainv - np.dot(Ainv, np.dot(Umat, np.dot(middle_part, np.dot(Vmat, Ainv))))

            
        L = np.linalg.cholesky(wmat)

        the_model = parsed["mu_bins"][stan_data["z_inds"][i]] + np.dot(stan_data["d_mag_d_filt"][i], parsed["dZPs"]) + np.dot(stan_data["dmagdcoeff"][i], parsed["coeff"])
        #the_model = parsed["mu_bins"][stan_data["z_inds"][i]] + parsed["dZPs"][stan_data["filt_inds"][i]] + np.dot(stan_data["dmagdcoeff"][i], parsed["coeff"])
                                          
        resid = stan_data["mags"][i] - the_model
        resid_norm.extend(np.dot(resid, L))

    priors = parsed["dZPs"]/0.005
    return np.concatenate((np.array(resid_norm), priors))


def run_fit(fit_coeff, fitdZP, outputsuffix):
    P, F, Cmat_no_model = miniLM_new(ministart = np.zeros(stan_data["NCoeff"] + stan_data["Nz"] + stan_data["NFilt"], dtype=np.float64),
                                     miniscale = unparseP(dict(coeff = [fit_coeff]*stan_data["NCoeff"], mu_bins = [1.]*stan_data["Nz"], dZPs = [fitdZP]*stan_data["NFilt"])),
                                     residfn = residfn,
                                     passdata = stan_data,
                                     verbose = True, maxiter = 1)
    
    mu_mat_no_model = Cmat_no_model[:stan_data["Nz"], :stan_data["Nz"]]
    comb_mat = np.zeros([len(mu_mat) + 1, len(mu_mat)], dtype=np.float64)
    comb_mat[0] = other_data["milliz_list"]
    comb_mat[0] /= 1000.
    
    comb_mat[1:] = mu_mat_no_model
    
    save_img(comb_mat, "comb_mat" + outputsuffix + ".fits")
    


direct_inverse = False # Invert cmat rather than use Woodbury matrix identity


parser = argparse.ArgumentParser()

parser.add_argument("pickle")
parser.add_argument("--prism_bins", help="Number of prism bins in wavelength", default = 25, type=int)
parser.add_argument("--SNRMax", help="Use >10 & >5 & >5 for LC selection, rather than total S/N", default = 0, type=int)
parser.add_argument("--model_res", help="Model wavelength resolution", default = 9, type=int)
parser.add_argument("--gray_disp", help="Gray dispersion", default = 0.1, type=float)
parser.add_argument("--twins", help="Assume twins S/N scaling", default = 0, type=int)
parser.add_argument("--color_scatter_opt", help="Color scatter optical", default = 0.03, type=float)
parser.add_argument("--color_scatter_nir", help="Color scatter nir", default = 0.03, type=float)
parser.add_argument("--train", help="Include Model Training", default = 1, type=int)
parser.add_argument("--calib", help="Include Calibration", default = 1, type=int)


opts = parser.parse_args()

prism_bins = opts.prism_bins

if opts.twins == 1:
    opts.gray_disp = -1
    opts.color_scatter_opt = 0.001
    opts.color_scatter_nir = 0.001

    
suffix = "_" + ( "notrain"*(opts.train == 0) + "nocalib"*(opts.calib == 0) ) + "SNRMax=%i_res=%02i_bins=%03i_disp=%.3f_scatopt=%.3f_scatnir=%.3f" % (opts.SNRMax, opts.model_res, opts.prism_bins, opts.gray_disp, opts.color_scatter_opt, opts.color_scatter_nir)



SN_data = load_data()
stan_data, other_data = get_stan_data(SN_data)




P, F, Cmat = miniLM_new(ministart = np.zeros(stan_data["NCoeff"] + stan_data["Nz"] + stan_data["NFilt"], dtype=np.float64),
                        miniscale = unparseP(dict(coeff = [opts.train]*stan_data["NCoeff"], mu_bins = [1.]*stan_data["Nz"], dZPs = [opts.calib]*stan_data["NFilt"])),
                        residfn = residfn,
                        passdata = stan_data,
                        verbose = True, maxiter = 1)

parsed = parseP(P, stan_data)

if opts.train and opts.calib:
    parsed_uncs = parseP(np.sqrt(np.diag(Cmat)), stan_data)
    
    print("parsed", parsed)
    print("parsed_uncs", parsed_uncs)


    rest_2D = rest_model_1D_to_2D(parsed["coeff"], other_data = other_data)
    rest_2D_mag_uncs = rest_model_1D_to_2D(parsed_uncs["coeff"], other_data = other_data)


    rest_mag_uncs_max = rest_2D_mag_uncs[other_data["phase0_ind"], :]
    
    assert len(rest_mag_uncs_max) == len(other_data["model_rest"])
    
    for i in range(len(other_data["model_rest"])):
        print("model_unc ", ("%05i" % int(other_data["model_rest"][i])), rest_mag_uncs_max[i])


    for i in range(len(other_data["filt_names"])):
        print("ZP_unc ", other_data["filt_names"][i], parsed_uncs["dZPs"][i])
    


    save_img(rest_2D, "rest_2D.fits")
    rest_2D_flux = 10.**(-0.4*rest_2D)
    save_img(rest_2D_flux, "rest_2D_flux.fits")

    save_img(rest_2D_mag_uncs, "rest_2D_mag_uncs.fits")


mu_mat = Cmat[:stan_data["Nz"], :stan_data["Nz"]]

print("mu_mat", mu_mat.shape)

comb_mat = np.zeros([len(mu_mat) + 1, len(mu_mat)], dtype=np.float64)
comb_mat[0] = other_data["milliz_list"]
comb_mat[0] /= 1000.

comb_mat[1:] = mu_mat

save_img(comb_mat, "comb_mat_" + suffix + ".fits")


#run_fit(fit_coeff = 0, fitdZP = 1, outputsuffix = suffix + "_no_model")
#run_fit(fit_coeff = 0, fitdZP = 0, outputsuffix = suffix + "_stat_only")
