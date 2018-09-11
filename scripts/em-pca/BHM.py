import pystan
import numpy as np
import cPickle as pickle

def make_sim_data():
    nsne = 1000

    

nbins = 50
[parameters, data, settings] = pickle.load(open("results.pickle", 'rb'))

for key in parameters:
    print "parameters ", key

for key in data:
    print "data ", key

print len(parameters["LC_fit_Cmat"])
print parameters["LC_fit_Cmat"][0].shape
print len(parameters["est_proj"])
print len(parameters["est_proj"][0])

stan_data = dict(nsne = 0, nbins = nbins, npar = 15, obs_par = [], obs_Cmat = [], zbins = [])

bin_boundaries = np.exp(np.linspace(np.log(np.amin(data["redshifts"]) - 1e-6 + 1.),
                                    np.log(np.amax(data["redshifts"]) + 1e-6 + 1.),
                                    nbins + 1)) - 1.

print bin_boundaries

for i in range(len(parameters["est_proj"])):
    if np.random.random() < 0.1:
        LC_par_Cmat = parameters["LC_fit_Cmat"][i][1:,1:]

        np.linalg.inv(LC_par_Cmat)
        np.linalg.cholesky(LC_par_Cmat)

        stan_data["nsne"] += 1
        stan_data["obs_par"].append(parameters["est_proj"][i])
        stan_data["obs_Cmat"].append(LC_par_Cmat)

        inds = np.where(data["redshifts"][i] >= bin_boundaries)
        gt_bins = bin_boundaries[inds]
        this_bin = len(gt_bins) - 1
        assert this_bin > -1 and this_bin < nbins

        stan_data["zbins"].append(this_bin)

print stan_data

stan_model = """
data {
    int nsne;
    int nbins;
    int npred; // Number of predictors, so 15 projections - 1 (magnitude) + 1 (color) = 15

    int zbins [nsne]; // Starts at 0!

    vector [npred + 1] obs_par [nsne];
    cov_matrix [npred + 1] obs_Cmat [nsne];
}

parameters {
    vector [npred] true_pred [nsne];
    real <lower = 0.02> sig_int;

    vector <lower = 0.25> [npred] R_npred;
    vector [npred] npred_star [nbins];
    vector [npred] coeffs;

    vector [nbins] mu;
}

model {
    matrix [npred + 1, npred + 1] obs_Cmat_with_int [nsne];
    vector [npred + 1] model_par [nsne];

    for (i in 1:nsne) {
        obs_Cmat_with_int[i] <- obs_Cmat[i];
        obs_Cmat_with_int[i][1,1] <- obs_Cmat[i][1,1] + sig_int^2;

        for (j in 1:npred) {
            model_par[i][j + 1] = true_par[i][j]; 
        }
        model_par[i][1] = mu + coeffs * true_pred[i];

        obs_par[i] ~ multi_normal(model_par[i], obs_Cmat_with_int[i]);
        true_par[i] ~ normal(par_star[zbins[i] + 1], R_par);
    }

    sig_int ~ normal(0, 0.1);
}
"""

sm = pystan.StanModel(model_code=stan_model)
fit = sm.sampling(data=stan_data, iter=1000, chains=4, refresh = 10)

print fit
