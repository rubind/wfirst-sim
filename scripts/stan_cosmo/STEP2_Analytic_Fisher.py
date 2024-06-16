import multiprocessing
multiprocessing.set_start_method("fork")
import pystan
import numpy as np
import pickle
import sys
import matplotlib.pyplot as plt
import extinction
from DavidsNM import save_img, miniLM_new
from scipy.interpolate import RectBivariateSpline
import tqdm


def load_data():
    SN_data = pickle.load(open(sys.argv[1], 'rb'))

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
    
def get_stan_data(SN_data):
    other_data = dict(filt_names = [], lambs_dict = dict(R062 = 6200., Z087 = 8700., Y106 = 10600., J129 = 12900., H158 = 15800., F184 = 18400., K213 = 21300., g = 4800., r = 6200., i = 7700., z = 8700.),
                      milliz_list = [], model_phase = np.linspace(-15., 45., 10), model_rest = np.exp(np.linspace(np.log(3000.), np.log(16000.), 9)))

    other_data["phase0_ind"] = np.argmin(np.abs(other_data["model_phase"]))
    other_data["restB_ind"] = np.argmin(np.abs(other_data["model_rest"] - 4400.))
    other_data["restV_ind"] = np.argmin(np.abs(other_data["model_rest"] - 5500.))

    
    stan_data = dict(NSNe = 0,
                     redshifts = [],
                     daymaxes = [],
                     rest_lambs = [],
                     true_phases = [],
                     color_law = [],
                     mags = [],
                     dmags = [],
                     z_inds = [],
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
                   
        if NPts > 10 and np.sqrt(sum(SNRs**2.)) > 40:
            stan_data["mags"].append(-2.5*np.log10(SN_LC["true_fluxes"][good_inds]))
            stan_data["dmags"].append((2.5/np.log(10.)) * np.abs(SN_LC["dfluxes"][good_inds]/SN_LC["true_fluxes"][good_inds]))

            this_daymax = SN_data["SN_table"]["daymaxes"][SN_ind]
            
            stan_data["redshifts"].append(this_z)
            stan_data["daymaxes"].append(this_daymax)
            

            for filt in np.unique(SN_LC["filts"][good_inds]):
                if other_data["filt_names"].count(filt) == 0:
                    other_data["filt_names"].append(filt)

            milliz = int(np.around(1000.*this_z))
            if other_data["milliz_list"].count(milliz) == 0:
                other_data["milliz_list"].append(milliz)
            
            stan_data["z_inds"].append(other_data["milliz_list"].index(milliz))        
            stan_data["filt_inds"].append([other_data["filt_names"].index(item) for item in SN_LC["filts"][good_inds]])
            
            stan_data["rest_lambs"].append(this_rest_lambs)
            stan_data["true_phases"].append((SN_LC["dates"][good_inds] - this_daymax)/(1. + this_z))
            stan_data["color_law"].append(extinction.ccm89(this_rest_lambs, a_v = 3.1, r_v = 3.1))
            
            stan_data["NSNe"] += 1

            
    stan_data["NFilt"] = len(other_data["filt_names"])
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

    
    print("stan_data", stan_data)
    print("other_data", other_data)
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
    

def residfn(P, wrapped_data):
    stan_data = wrapped_data[0]
    parsed = parseP(P, stan_data)
    
    resid_norm = []
    
    for i in tqdm.trange(stan_data["NSNe"]):
        cmat = np.diag(stan_data["dmags"][i]**2.)
        cmat += 0.1**2.
        cmat += np.outer(stan_data["color_law"][i]*0.5, stan_data["color_law"][i]*0.5) # +- 0.5 mag E(B-V)
        
        for filt in np.unique(stan_data["filt_inds"][i]):
            filt_mask = (stan_data["filt_inds"][i] == filt)
            cmat += np.outer(filt_mask*0.03, filt_mask*0.03)

        wmat = np.linalg.inv(cmat)
        L = np.linalg.cholesky(wmat)

        the_model = parsed["mu_bins"][stan_data["z_inds"][i]] + parsed["dZPs"][stan_data["filt_inds"][i]] + np.dot(stan_data["dmagdcoeff"][i], parsed["coeff"])
        resid = stan_data["mags"][i] - the_model
        resid_norm.extend(np.dot(resid, L))

    priors = parsed["dZPs"]/0.005
    return np.concatenate((np.array(resid_norm), priors))





SN_data = load_data()
stan_data, other_data = get_stan_data(SN_data)


P, F, Cmat = miniLM_new(ministart = np.zeros(stan_data["NCoeff"] + stan_data["Nz"] + stan_data["NFilt"], dtype=np.float64),
                        miniscale = np.ones(stan_data["NCoeff"] + stan_data["Nz"] + stan_data["NFilt"], dtype=np.float64),
                        residfn = residfn,
                        passdata = stan_data,
                        verbose = True, maxiter = 1)

parsed = parseP(P, stan_data)
parsed_uncs = parseP(np.sqrt(np.diag(Cmat)), stan_data)

print("parsed", parsed)
print("parsed_uncs", parsed_uncs)


rest_2D = rest_model_1D_to_2D(parsed["coeff"], other_data = other_data)
rest_2D_mag_uncs = rest_model_1D_to_2D(parsed_uncs["coeff"], other_data = other_data)

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

save_img(comb_mat, "comb_mat.fits")
