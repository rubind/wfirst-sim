from numpy import *
import pystan
from DavidsNM import save_img
import sys
import os
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/synth_dataset/")
import synth_functions
from astropy import cosmology
from extinction import fm07
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import cPickle as pickle
import commands
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import time


def get_EVs_and_color_law(filt, z, filtfns, phases, cosmo, filtwaves, nEV, meanfn, compfns):
    lambs = arange(filtwaves[filt][0]*0.8, filtwaves[filt][1]*1.25, 10)
    ab_mags_mean = synth_functions.get_AB_mags(true_norm = 1., true_proj = zeros(nEV, dtype=float64), true_AV = 0, z = z, filtfn = filtfns[filt],
                                               phases = phases, cosmo = cosmo, lambs = lambs, ncomp = nEV,
                                               meanfn = meanfn, compfns = compfns)*1e14

    color_law = -2.5*log10(synth_functions.get_AB_mags(true_norm = 1., true_proj = zeros(nEV, dtype=float64), true_AV = 1, z = z, filtfn = filtfns[filt],
                                                       phases = phases, cosmo = cosmo, lambs = lambs, ncomp = nEV,
                                                       meanfn = meanfn, compfns = compfns)*1e14/ab_mags_mean)
    
    tmp_eval = []
    for i in range(nEV):
        tmp_proj = zeros(nEV, dtype=float64)
        tmp_proj[i] = 1

        ab_mags = synth_functions.get_AB_mags(true_norm = 1., true_proj = tmp_proj, true_AV = 0, z = z, filtfn = filtfns[filt],
                                              phases = phases, cosmo = cosmo, lambs = lambs, ncomp = nEV,
                                              meanfn = meanfn, compfns = compfns)*1e14
        ab_mags -= ab_mags_mean

        tmp_eval.append(ab_mags)
    return ab_mags_mean, array(tmp_eval), color_law


def get_EVs_and_color_law_prism(z, cosmo, nEV, meanfn, compfns):
    if z > 0.1:
        lambs = exp(
            arange(log(max(6000., 3300.*(1. + z))),
                   log(min(18000., 8600.*(1. + z))), 0.01))
    else:
        lambs = exp(arange(log(3300*(1. + z)), log(8600*(1. + z)), 0.01))



    ab_spec_mean = synth_functions.get_spectrum(true_norm = 1., true_proj = zeros(nEV, dtype=float64), true_AV = 0, restlambs = lambs/(1. + z), ncomp = nEV, meanfn = meanfn, compfns = compfns)


    color_law = -2.5*log10(synth_functions.get_spectrum(true_norm = 1., true_proj = zeros(nEV, dtype=float64), true_AV = 1., restlambs = lambs/(1. + z), ncomp = nEV, meanfn = meanfn, compfns = compfns)/
                           ab_spec_mean)
    
    tmp_eval = []
    for i in range(nEV):
        tmp_proj = zeros(nEV, dtype=float64)
        tmp_proj[i] = 1

        ab_spec = synth_functions.get_spectrum(true_norm = 1., true_proj = tmp_proj, true_AV = 0, restlambs = lambs/(1. + z), ncomp = nEV, meanfn = meanfn, compfns = compfns)
        ab_spec -= ab_spec_mean

        tmp_eval.append(ab_spec)
    return ab_spec_mean, array(tmp_eval), color_law


def make_sim_SN(i, SN_data, stan_data, filts, spec_SNR = None):
    stan_data["nsne"] += 1

    if i != None:
        phases = (SN_data["SN_observations"][i]["dates"] - SN_data["SN_table"]["daymaxes"][i])/(1. + SN_data["SN_table"]["redshifts"][i])
        stan_data["redshifts"].append(SN_data["SN_table"]["redshifts"][i])
    else:
        phases = arange(-10, 31., 4.)
        stan_data["redshifts"].append(0.05)

    stan_data["lowind"].append(stan_data["nobs"])
    stan_data["highind"].append(stan_data["nobs"] - 1)


    true_AV = random.exponential()*0.3
    true_dM = random.normal()*args.sigint
    true_ind = random.randint(nx0)
    true_proj = projs[true_ind]
    true_norm = smoothed_norms[true_ind]

    stan_data["true_AVs"].append(true_AV)
    stan_data["true_dMs"].append(true_dM)
    stan_data["true_inds"].append(true_ind)
    stan_data["true_projs"].append(true_proj)
    stan_data["true_norms"].append(true_norm)


    for filt in filts:
        if i != None:
            inds = where((SN_data["SN_observations"][i]["filts"] == filt)*filts_are_okay_mask*(phases > -10)*(phases < 30))
            these_SNRs = clip(SNRs[inds], 0.1, 10000)
            print filt, len(inds[0]), these_SNRs

            stan_data["nobs"] += len(inds[0])
            stan_data["highind"][-1] += len(inds[0])


            the_mean, the_EVs, the_color_law = get_EVs_and_color_law(filt, stan_data["redshifts"][-1], filtfns, phases[inds], cosmo, filtwaves, nEV, meanfn, compfns)

        else:
            these_SNRs = ones(len(phases), dtype=float64)*25.

            stan_data["nobs"] += len(phases)
            stan_data["highind"][-1] += len(phases)

            the_mean, the_EVs, the_color_law = get_EVs_and_color_law(filt, stan_data["redshifts"][-1], filtfns, phases, cosmo, filtwaves, nEV, meanfn, compfns)

        stan_data["mean_eval"].extend(the_mean)
        stan_data["color_law"].extend(the_color_law)
        for j in range(nEV):
            stan_data["EVs"][j].extend(the_EVs[j])

        true_flux = true_norm*(the_mean + dot(true_proj, the_EVs))*10.**(-0.4*(true_dM + true_AV*the_color_law))
        true_err = abs(true_flux/these_SNRs)

        stan_data["obs_fluxes"].extend(true_flux + random.normal(size = len(true_flux))*true_err)
        stan_data["obs_dfluxes"].extend(true_err)

    if spec_SNR != None:
        assert i == None

        the_mean, the_EVs, the_color_law = get_EVs_and_color_law_prism(stan_data["redshifts"][-1],
                                                                       cosmo, nEV, meanfn, compfns)
        stan_data["nobs"] += len(the_mean)
        stan_data["highind"][-1] += len(the_mean)

        stan_data["mean_eval"].extend(the_mean)
        stan_data["color_law"].extend(the_color_law)
        for j in range(nEV):
            stan_data["EVs"][j].extend(the_EVs[j])

        true_flux = true_norm*(the_mean + dot(true_proj, the_EVs))*10.**(-0.4*(true_dM + true_AV*the_color_law))
        true_err = abs(true_flux/spec_SNR)

        stan_data["obs_fluxes"].extend(true_flux + random.normal(size = len(true_flux))*true_err)
        stan_data["obs_dfluxes"].extend(true_err)

    return stan_data


stan_model = """
data {
    int nsne;
    int nEV;
    int nobs; // All observations, imaging+spectra, all SNe
    int nredcoeff; // E.g., 2 for linear variation

    int lowind [nsne];
    int highind [nsne];
    matrix [nobs, nEV] EVs;
    vector [nredcoeff] redcoeff [nsne];
    vector [nobs] mean_eval;
    vector [nobs] color_law;
    
    vector [nobs] obs_fluxes;
    vector [nobs] obs_dfluxes;
}

parameters {
    matrix [nEV, nredcoeff] xstar;
    matrix <lower = 0> [nEV, nredcoeff] Rx;
    vector [nEV] true_x [nsne];
    vector <lower = 0> [nsne] true_AV;
    vector [nsne] dM;
    row_vector <lower = 0> [nredcoeff] tau_AV;
}


model {
    vector [nobs] mod_fluxes;
    vector [nEV] xstar_by_SN [nsne];
    vector [nEV] Rx_by_SN [nsne];
    vector [nsne] tau_AV_by_SN;

    for (i in 1:nsne) {
        xstar_by_SN[i] <- xstar * redcoeff[i];
        Rx_by_SN[i] <- Rx * redcoeff[i];
        tau_AV_by_SN[i] <- tau_AV * redcoeff[i];
    }

    for (i in 1:nsne) {
        mod_fluxes[lowind[i]+1:highind[i]+1] <- (EVs[lowind[i]+1:highind[i]+1]*(true_x[i] + xstar_by_SN[i]) + mean_eval[lowind[i]+1:highind[i]+1]) .*
                                                exp(-0.4*log(10.)*(   color_law[lowind[i]+1:highind[i]+1]*true_AV[i] + dM[i]   ));
        true_x[i] ~ normal(0, Rx_by_SN[i]);
        true_AV[i] ~ exponential(1./(tau_AV_by_SN[i]));
    }

    obs_fluxes ~ normal(mod_fluxes, obs_dfluxes);
    
    for (i in 1:nredcoeff) {
        tau_AV[i] ~ normal(0.3, 0.1);
        Rx[i] ~ normal(0, 10);
        xstar[i] ~ normal(0, 10);
    }
}
"""


def initfn():
    return dict(xstar = random.normal(size = [nEV, 2])*0.2,
                Rx = random.random(size = [nEV, 2]) + 0.5,
                true_x = stan_data["true_projs"] - median(stan_data["true_projs"], axis = 0) + random.normal(size = [stan_data["nsne"], nEV])*0.05,
                true_AV = stan_data["true_AVs"],
                true_dM = random.normal(size = stan_data["nsne"])*0.01,
                tau_AV = 0.3 + random.random(size = 2)*0.02)



meanfn, norms, compfns, projs, smoothed_norms, NA = synth_functions.get_data()
nEV = len(compfns)
nx0 = len(smoothed_norms)
filtfns, filtwaves, filtnames = synth_functions.get_filts()
print "projs.shape", projs.shape
print "filtwaves", filtwaves
print "nEV", nEV
print "nx0", nx0

cosmo = cosmology.FlatLambdaCDM(Om0 = 0.3, H0 = 70.)


parser = argparse.ArgumentParser()
parser.add_argument("--sigint", type = float, help = "intrinsic dispersion", default = 0.08)
parser.add_argument("--fast", type = int, help = "fast (for testing)", default = 0)
parser.add_argument("--pickle", type = str, help = "pickle file")
args = parser.parse_args()

print time.asctime()
SN_data = pickle.load(open(args.pickle, 'rb'))
print time.asctime()



stan_data = dict(nsne = 0, nredcoeff = 2, nobs = 0, nEV = 14,
                 redshifts = [], redshifts_has_IFC = [], has_IFC = [],
                 lowind = [], highind = [],
                 EVs = [[] for i in range(nEV)], mean_eval = [], color_law = [],

                 true_AVs = [], true_dMs = [], true_inds = [], true_projs = [], true_norms = [],

                 obs_fluxes = [], obs_dfluxes = [])


okay_filters = {}
for z in unique(SN_data["SN_table"]["redshifts"]):
    okay_filters[z] = [item for item in filtwaves if filtwaves[item][0]/(1. + z) > 3300. and filtwaves[item][1]/(1. + z) < 8600.]

print "okay_filters ", okay_filters

for ii in range(800 - args.fast * 700):
    print ii
    stan_data = make_sim_SN(None, SN_data, stan_data, filts = ["g", "r", "i"], spec_SNR = 20.)

pdf = PdfPages("SNe_and_selection.pdf")

if args.fast:
    SN_data["nsne"] = 200

for i in range(SN_data["nsne"]):
    this_z = SN_data["SN_table"]["redshifts"][i]

    if len(SN_data["SN_observations"][i]["fluxes"]) > 4:
        has_IFC = len(SN_data["SN_observations"][i]["IFS_dates"]) > 0
        
        SNRs = array(SN_data["SN_observations"][i]["fluxes"])/array(SN_data["SN_observations"][i]["dfluxes"])
        phases = (SN_data["SN_observations"][i]["dates"] - SN_data["SN_table"]["daymaxes"][i])/(1. + this_z)


        filts_are_okay_mask = [okay_filters[this_z].count(filtitem) for filtitem in SN_data["SN_observations"][i]["filts"]]
        inds = where((SNRs > 8)*filts_are_okay_mask*(phases > -10)*(phases < 30))

        filts = unique(SN_data["SN_observations"][i]["filts"][inds])
        if ((len(filts) > 2) and (sum(SNRs > 8) > 4)):# or has_IFC:
            print "Using SN ", i, " of ", SN_data["nsne"]
            selected = 1

            stan_data = make_sim_SN(i, SN_data, stan_data, filts)

            stan_data["has_IFC"].append(int(has_IFC))
            if has_IFC:
                stan_data["redshifts_has_IFC"].append(SN_data["SN_table"]["redshifts"][i])

        else:
            print "Not using SN ", i
            selected = 0

        if random.random() < 0.005:
            plt.figure(figsize = (7,5))
            for filt in unique(SN_data["SN_observations"][i]["filts"]):
                inds = where(SN_data["SN_observations"][i]["filts"] == filt)
                plt.errorbar(SN_data["SN_observations"][i]["dates"][inds], SN_data["SN_observations"][i]["fluxes"][inds], yerr = SN_data["SN_observations"][i]["dfluxes"][inds], fmt = '.', label = filt)

            for IFCdate in SN_data["SN_observations"][i]["IFS_dates"]:
                plt.axvline(IFCdate, color = 'k')

            plt.legend(loc = 'best')
            plt.title((1 - selected)*"NOT " + "Selected, z=%.3f" % (SN_data["SN_table"]["redshifts"][i]))
            plt.axhline(0, color = 'gray')
            

            pdf.savefig(plt.gcf())
            plt.close()
pdf.close()

print "nobs ", stan_data["nobs"]
print "nsne ", stan_data["nsne"]

print stan_data["lowind"][:10]
print stan_data["highind"][:10]


stan_data["mean_eval"] = array(stan_data["mean_eval"])
stan_data["color_law"] = array(stan_data["color_law"])
stan_data["EVs"] = transpose(array(stan_data["EVs"]))
stan_data["obs_fluxes"] = array(stan_data["obs_fluxes"])
stan_data["obs_dfluxes"] = array(stan_data["obs_dfluxes"])

for key in ["color_law", "mean_eval", "EVs", "obs_fluxes"]:
    inds = where(isnan(stan_data[key]) + isinf(stan_data[key]))
    stan_data[key][inds] = 0

inds = where(isnan(stan_data["obs_dfluxes"]) + isinf(stan_data["obs_dfluxes"]))
stan_data["obs_dfluxes"][inds] = stan_data["obs_fluxes"].max()


save_img(stan_data["EVs"], "EVs.fits")
save_img(stan_data["obs_fluxes"], "obs_fluxes.fits")
save_img(stan_data["obs_dfluxes"], "obs_dfluxes.fits")
save_img(stan_data["color_law"], "color_law.fits")

for key in ["mean_eval", "color_law", "EVs"]:
    print key, stan_data[key].shape

plt.hist(stan_data["redshifts"], bins = arange(0., 2.55, 0.05), label = str(len(stan_data["redshifts"])))
plt.hist(stan_data["redshifts_has_IFC"], bins = arange(0., 2.55, 0.05), label = str(len(stan_data["redshifts_has_IFC"])))
plt.legend(loc = 'best')
plt.savefig("redshifts_of_SNe_selected.pdf")
plt.close()


stan_data["redshifts"] = array(stan_data["redshifts"])
stan_data["redcoeff"] = transpose(array([
    1 - 1./(1. + stan_data["redshifts"]),
    1./(1. + stan_data["redshifts"])]))
    
print "redcoeff", stan_data["redcoeff"].shape

pickle.dump(stan_data, open("input.pickle", 'wb'))

sm = pystan.StanModel(model_code=stan_model)

fit = sm.sampling(data=stan_data, iter=1000, chains=4, refresh = 1, init = initfn)
fit_params = fit.extract(permuted = True)

print fit

pickle.dump([stan_data, fit_params], open("fit_results_and_input.pickle", 'wb'))

#print commands.getoutput("~/Dropbox/imessage.sh davidrubinlbl@gmail.com 'Done running!'")
