from numpy import *
import pystan
from DavidsNM import save_img
import sys
sys.path.append("../synth_dataset/")
import synth_functions
from astropy import cosmology
from extinction import fm07
import matplotlib.pyplot as plt
import cPickle as pickle
import commands
import argparse


stan_model = """
data {
    int nsne;
    int nEV;
    int ndim;
    int nx0;
    int nspecdim;
    int nspectra;
    int nredshiftbin;
    int redshiftbin [nsne];

    matrix [nEV, ndim] EVs [nsne];
    vector [ndim] mean_eval [nsne];
    vector [nspecdim] dspec_dAV;
    
    vector [ndim] color_law [nsne];
    vector [ndim] obs_fluxes [nsne];
    vector [ndim] obs_dfluxes [nsne];

    int spec_SN_inds [nspectra];

    vector [nspecdim] obs_spectra [nspectra];
    vector [nspecdim] obs_dspectra [nspectra];
    vector [nspecdim] mean_spec_eval;
    matrix [nEV, nspecdim] comp_spec_evals;
}

parameters {
    row_vector [nEV] xstar [nredshiftbin];
    vector <lower = 0> [nEV] Rx [nredshiftbin];
    row_vector [nEV] true_x [nsne];
    vector <lower = 0> [nsne] true_AV;
    vector [nsne] dM;
    real <lower = 0> tau_AV [nredshiftbin];
}


model {
    vector [ndim] mod_fluxes [nsne];

    for (i in 1:nsne) {

        mod_fluxes[i] <- ((true_x[i] + xstar[redshiftbin[i] + 1]) * EVs[i] + mean_eval[i]')' .* exp(-0.4*log(10.)*(   color_law[i]*true_AV[i] + dM[i]   ));
        obs_fluxes[i] ~ normal(mod_fluxes[i], obs_dfluxes[i]);
        true_x[i] ~ normal(0, Rx[redshiftbin[i] + 1]);
        true_AV[i] ~ exponential(1./tau_AV[redshiftbin[i] + 1]);
    }

    for (i in 1:nspectra) {
        obs_spectra[i] ~ normal(
                               ((true_x[spec_SN_inds[i] + 1] + xstar[redshiftbin[i] + 1]) * comp_spec_evals + mean_spec_eval' )' .* 
                               exp(-0.4*log(10.) * (dspec_dAV * true_AV[spec_SN_inds[i] + 1] + dM[spec_SN_inds[i] + 1]))
                               , obs_dspectra[i]);
    }
    
    for (i in 1:nredshiftbin) {
        tau_AV[i] ~ normal(0.3, 0.1);
        Rx[i] ~ normal(0, 10);
        xstar[i] ~ normal(0, 10);
    }
}
"""


def initfn():
    return dict(xstar = random.normal(size = nEV)*0.2,
                Rx = random.random(size = nEV) + 0.5,
                true_x = true_projs + random.normal(size = true_projs.shape)*0.05,
                true_AV = true_AVs,
                true_dM = random.normal(size = nsne)*0.01,
                M = random.normal(),
                tau_AV = 0.3 + random.random()*0.02)



meanfn, norms, compfns, projs, smoothed_norms, NA = synth_functions.get_data()
filtfns, filtwaves, filtnames = synth_functions.get_filts()
print "projs.shape", projs.shape


z = 1.0
cosmo = cosmology.FlatLambdaCDM(Om0 = 0.3, H0 = 70.)

dates = arange(-15., 70., 5.)
print dates, dates/(1. + z), len(dates)

filts_to_run = [item for item in filtwaves if filtwaves[item][0]/(1. + z) > 3300. and filtwaves[item][1]/(1. + z) < 8600.]

print "filts_to_run", filts_to_run


parser = argparse.ArgumentParser()
parser.add_argument("--specfrac", type = float, help = "fraction of SNe with spectrscopy")
parser.add_argument("--maxSNR", type = float, help = "imaging SNR at peak")
parser.add_argument("--maxspecSNR", type = float, help = "spectroscopic max SNR per 3000 km/s")
parser.add_argument("--nsne", type = int, help = "number of SNe")
parser.add_argument("--sigint", type = float, help = "intrinsic dispersion")
args = parser.parse_args()


ndim = len(filts_to_run)*len(dates) # Number of datapoints
nEV = len(projs[0]) # Number of EVs
nsne = args.nsne # Number of SNe
nx0 = len(projs) # Number of archetypes
maxSNR = args.maxSNR
spec_frac = args.specfrac
maxspecSNR = args.maxspecSNR
sigint = args.sigint

spec_lambs = exp(arange(log(3300.), log(8600.), 0.01))
nspecdim = len(spec_lambs)

stan_data = dict(nsne = nsne, nEV = nEV, ndim = ndim, nx0 = nx0)
comp_evals = [] # Flux due to one unit of each comp

comp_spec_evals = []

for i in range(nEV):
    tmp_eval = []
    tmp_proj = zeros(nEV, dtype=float64)
    tmp_proj[i] = 1
    mean_eval = []

    mean_spec_eval = synth_functions.get_spectrum(true_norm = 1., true_proj = zeros(nEV, dtype=float64), true_AV = 0, restlambs = spec_lambs, ncomp = nEV, meanfn = meanfn, compfns = compfns)*spec_lambs
    the_spec = synth_functions.get_spectrum(true_norm = 1., true_proj = tmp_proj, true_AV = 0, restlambs = spec_lambs, ncomp = nEV, meanfn = meanfn, compfns = compfns)*spec_lambs
    comp_spec_evals.append(the_spec - mean_spec_eval)

    for filt in filts_to_run:
        lambs = arange(filtwaves[filt][0]*0.8, filtwaves[filt][1]*1.25, 10)
        ab_mags = synth_functions.get_AB_mags(true_norm = 1., true_proj = tmp_proj, true_AV = 0, z = 1.0, filtfn = filtfns[filt],
                                              phases = dates/(1. + z), cosmo = cosmo, lambs = lambs, ncomp = nEV,
                                              meanfn = meanfn, compfns = compfns)*1e14
        ab_mags_mean = synth_functions.get_AB_mags(true_norm = 1., true_proj = zeros(nEV, dtype=float64), true_AV = 0, z = 1.0, filtfn = filtfns[filt],
                                                   phases = dates/(1. + z), cosmo = cosmo, lambs = lambs, ncomp = nEV,
                                                   meanfn = meanfn, compfns = compfns)*1e14
        ab_mags -= ab_mags_mean

        tmp_eval.extend(ab_mags)
        mean_eval.extend(ab_mags_mean)

    comp_evals.append(tmp_eval)
comp_evals = array(comp_evals)
mean_eval = array(mean_eval)
comp_spec_evals = array(comp_spec_evals)


print "comp_evals", comp_evals.shape

save_img(comp_evals, "comp_evals.fits")
save_img(mean_eval, "mean_eval.fits")
save_img(comp_spec_evals, "comp_spec_evals.fits")


dspec_dAV = -2.5*log10(
    synth_functions.get_spectrum(true_norm = 1., true_proj = median(projs, axis = 0), true_AV = 1, restlambs = spec_lambs, ncomp = nEV, meanfn = meanfn, compfns = compfns)/
    synth_functions.get_spectrum(true_norm = 1., true_proj = median(projs, axis = 0), true_AV = 0, restlambs = spec_lambs, ncomp = nEV, meanfn = meanfn, compfns = compfns))

dmag_dAV = []
for filt in filts_to_run:
    lambs = arange(filtwaves[filt][0]*0.8, filtwaves[filt][1]*1.25, 10)
    ab_mags1 = synth_functions.get_AB_mags(true_norm = 1, true_proj = median(projs, axis = 0), true_AV = 1, z = 1.0, filtfn = filtfns[filt],
                                           phases = dates/(1. + z), cosmo = cosmo, lambs = lambs, ncomp = nEV, compfns = compfns, meanfn = meanfn)*1e14
    ab_mags0 = synth_functions.get_AB_mags(true_norm = 1, true_proj = median(projs, axis = 0), true_AV = 0, z = 1.0, filtfn = filtfns[filt],
                                           phases = dates/(1. + z), cosmo = cosmo, lambs = lambs, ncomp = nEV, compfns = compfns, meanfn = meanfn)*1e14
    dmag_dAV.extend(-2.5*log10(ab_mags1/ab_mags0))

dmag_dAV = array(dmag_dAV)
plt.plot(dmag_dAV)
plt.savefig("dmag_dAV.pdf")
plt.close()

plt.plot(dspec_dAV)
plt.savefig("dspec_dAV.pdf")
plt.close()


obs_fluxes = []
obs_dfluxes = []
true_projs = []

true_AVs = random.exponential(size = nsne)*0.3
true_dMs = random.normal(size = nsne)*sigint

obs_spectra = []
obs_dspectra = []
spec_SN_inds = []
nspectra = 0

for i in range(nsne):
    print "*"
    this_ind = random.randint(nx0)
    this_proj = projs[this_ind]
    true_projs.append(this_proj)


    this_flux = []
    for filt in filts_to_run:
        lambs = arange(filtwaves[filt][0]*0.8, filtwaves[filt][1]*1.25, 10)
        ab_mags = synth_functions.get_AB_mags(true_norm = smoothed_norms[this_ind]*10.**(-0.4*true_dMs[i]),
                                              true_proj = this_proj, true_AV = true_AVs[i], z = 1.0, filtfn = filtfns[filt],
                                              phases = dates/(1. + z), cosmo = cosmo, lambs = lambs, ncomp = nEV,
                                              meanfn = meanfn, compfns = compfns)*1e14
        this_flux.extend(ab_mags)
    
    this_flux = array(this_flux)
    this_err = abs(this_flux.max()) / maxSNR
    obs_dfluxes.append(ones(ndim, dtype=float64)*this_err)
    obs_fluxes.append(this_flux + random.normal(size = ndim)*obs_dfluxes[-1])

    if random.random() < spec_frac or i == 0:
        print "Making spectrum!"
    
        this_spectrum = synth_functions.get_spectrum(true_norm = smoothed_norms[this_ind]*10.**(-0.4*true_dMs[i]),
                                                     true_proj = this_proj, true_AV = true_AVs[i], restlambs = spec_lambs,
                                                     ncomp = nEV, meanfn = meanfn, compfns = compfns)*spec_lambs
        this_err = abs(max(this_spectrum))/ maxspecSNR

        obs_spectra.append(this_spectrum + random.normal(size = nspecdim)*this_err)
        obs_dspectra.append(ones(nspecdim, dtype=float64)*this_err)
        spec_SN_inds.append(i)
        nspectra += 1



true_projs = array(true_projs)
obs_fluxes = array(obs_fluxes)
obs_dfluxes = array(obs_dfluxes)

obs_spectra = array(obs_spectra)
obs_dspectra = array(obs_dspectra)

save_img(obs_fluxes, "obs_fluxes.fits")
save_img(obs_spectra, "obs_spectra.fits")


stan_data["obs_fluxes"] = obs_fluxes
stan_data["obs_dfluxes"] = obs_dfluxes
stan_data["EVs"] = comp_evals
stan_data["mean_eval"] = mean_eval
stan_data["color_law"] = dmag_dAV
stan_data["projs_0"] = projs
stan_data["norms_0"] = smoothed_norms

stan_data["nspecdim"] = nspecdim
stan_data["obs_spectra"] = obs_spectra
stan_data["obs_dspectra"] = obs_dspectra
stan_data["spec_SN_inds"] = spec_SN_inds
stan_data["dspec_dAV"] = dspec_dAV
stan_data["mean_spec_eval"] = mean_spec_eval
stan_data["comp_spec_evals"] = comp_spec_evals
stan_data["nspectra"] = nspectra

sm = pystan.StanModel(model_code=stan_model)

fit = sm.sampling(data=stan_data, iter=1000, chains=4, refresh = 1, init = initfn)
fit_params = fit.extract(permuted = True)

print fit

pickle.dump([stan_data, fit_params], open("fit_results_NSNe=%i_maxSNR=%.1f_specfrac=%.3f_specSNR=%.1f_sigint=%.3f.pickle" % (nsne, maxSNR, spec_frac, maxspecSNR, sigint), 'wb'))

#print commands.getoutput("~/Dropbox/imessage.sh davidrubinlbl@gmail.com 'Done running!'")
