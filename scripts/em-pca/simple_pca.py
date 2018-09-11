import numpy as np
from scipy.interpolate import RectBivariateSpline
from DavidsNM import miniLM_new
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import time
import sys
import cPickle as pickle
import gzip
import multiprocessing as mp


def file_to_dict(flname):
    settings = {}

    f = open(flname, 'r')
    lines = f.read().split('\n')
    f.close()

    for line in lines:
        parsed = line.split("#")[0].split(None)
        if len(parsed) > 1:
            settings[parsed[0]] = eval(" ".join(parsed[1:]))
    return settings


def get_settings():
    settings = file_to_dict(sys.argv[1])

    data = {}
    parameters = {}

    try:
        print "Trying gzipped file"
        (survey, filter_info, EVs_with_filter, lambsteps, phases) = pickle.load(gzip.open(settings["survey_pickle"]))
    except:
        print "Nope!"
        (survey, filter_info, EVs_with_filter, lambsteps, phases) = pickle.load(open(settings["survey_pickle"]))


    print EVs_with_filter.shape
    for i in range(len(EVs_with_filter)):
        print i, "EVs_with_filter RMS", np.std(EVs_with_filter[i])

    settings["nsn"] = len(survey)
    settings["phases"] = phases
    settings["nph"] = len(settings["phases"])
    settings["nev"] = len(EVs_with_filter) - 1 # Doesn't include color law, hence the - 1

    settings["rest_nodes"] = np.exp(np.linspace(np.log(lambsteps.min()),
                                                np.log(lambsteps.max()),
                                                settings["nlm"])) # Filters have log width filt_spacing

    settings["EV_treatment"] = np.array(settings["EV_treatment"])

    parameters["est_EV_splines"] = []
    for i in range(settings["nev"]):
        if settings["EV_treatment"][i] == 2:
            parameters["est_EV_splines"].append(None)
        else:
            parameters["est_EV_splines"].append(RectBivariateSpline(lambsteps, settings["phases"], EVs_with_filter[i], kx = 3, ky = 3))

    
    parameters["color_law"] = RectBivariateSpline(lambsteps, settings["phases"], EVs_with_filter[-1], kx = 3, ky = 3)


    median_flux = np.median([np.median(survey[i]["lc"]["flux"]) for i in range(settings["nsn"])])
    print "median_flux ", median_flux

    data["redshifts"] = []
    data["dates"] = []
    data["rlambs"] = []
    data["fluxes"] = []
    data["invvars"] = []

    for i in range(settings["nsn"]):
        data["redshifts"].append(survey[i]["z"])
        
        unique_filts = list(np.sort(np.unique(survey[i]["lc"]["band"])))
        data["dates"].append(list(np.sort(np.unique(survey[i]["lc"]["time"]))))
        data["rlambs"].append([filter_info[item]/(1. + survey[i]["z"]) for item in unique_filts])
        
        
        assert len(survey[i]["lc"]["flux_orig"]) == len(data["dates"][-1])*len(data["rlambs"][-1])
        
        flux_grid = np.zeros([len(data["rlambs"][-1]), len(data["dates"][-1])], dtype=np.float64)
        dflux_grid = np.zeros([len(data["rlambs"][-1]), len(data["dates"][-1])], dtype=np.float64)
        
        for j in range(len(survey[i]["lc"]["flux_orig"])):
            ind_i, ind_j = unique_filts.index(survey[i]["lc"]["band"][j]), data["dates"][-1].index(survey[i]["lc"]["time"][j]),

            flux_grid[ind_i, ind_j] = survey[i]["lc"]["flux"][j]
            dflux_grid[ind_i, ind_j] = survey[i]["lc"]["flux_err"][j]

        data["fluxes"].append(flux_grid/median_flux)
        data["invvars"].append((dflux_grid/median_flux)**-2.)


    print "settings:"
    for key in settings:
        print key, settings[key]


    return parameters, data, settings

def get_initial_parameters(parameters, settings):

    parameters["est_EVs"] = np.random.normal(size = [sum(settings["EV_treatment"] == 2), settings["nlm"], settings["nph"]])
    parameters["est_proj"] = np.zeros([settings["nsn"], settings["nev"] + 1], dtype=np.float64)
    parameters["est_proj"][:, 0] = 0.

    parameters["est_daymax"] = np.zeros(settings["nsn"], dtype=np.float64)
    parameters["LC_fit_Cmat"] = [None for i in range(settings["nsn"])]

    parameters = get_splines(parameters, settings)

    return parameters

def get_splines(parameters, settings):
    ind = 0
    assert len(parameters["est_EVs"]) == sum(settings["EV_treatment"] == 2)

    for i in range(settings["nev"]):
        if settings["EV_treatment"][i] == 2:
            parameters["est_EV_splines"][i] = RectBivariateSpline(settings["rest_nodes"], settings["phases"], parameters["est_EVs"][ind], kx = 3, ky = 3)
            ind += 1

    return parameters


def modelfn(parameters, settings, sne_to_do = None):
    if sne_to_do == None:
        sne_to_do = range(settings["nsn"])
    else:
        try:
            sne_to_do[0] # list, array, whatever...
        except:
            sne_to_do = [sne_to_do]

    the_model = []

    for i in sne_to_do:
        z = data["redshifts"][i]
        phases = (data["dates"][i] - parameters["est_daymax"][i])/(1. + z)
        rlambs = data["rlambs"][i]

        this_model = parameters["est_EV_splines"][0](rlambs, phases, grid=True)
        for j in range(1, settings["nev"]): # Not including color
            this_model += parameters["est_proj"][i][j]*parameters["est_EV_splines"][j](rlambs, phases, grid=True)
        this_model *= 10.**(-0.4*(   parameters["est_proj"][i][0] + parameters["color_law"](rlambs, phases, grid=True) * parameters["est_proj"][i,-1]   ))

        the_model.append(this_model)
    return the_model
        

def E_pullfn(P, passdata):
    [i, parameters, data, settings] = passdata[0]

    parameters["est_daymax"][i] = P[0]
    parameters["est_proj"][i] = P[1:]

    the_model = modelfn(parameters, settings, sne_to_do = i)
    pulls = np.sqrt(data["invvars"][i])*(the_model[0] - data["fluxes"][i])
    
    return np.concatenate((pulls.flatten(),
                           P/100.))



def E_step(parameters, data, settings):
    print "Starting E step..."

    miniscale = 100*np.ones(2 + settings["nev"], dtype=np.float64)
    for i in range(settings["nev"]):
        if settings["EV_treatment"][i] == 0: # Ignore this eigenvector
            miniscale[i+1] = 0

    

    for i in range(settings["nsn"]):
        P, F, Cmat = miniLM_new(ministart = np.concatenate(([parameters["est_daymax"][i]], parameters["est_proj"][i])),
                                miniscale = miniscale,
                                residfn = E_pullfn, passdata = [i, parameters, data, settings], verbose = False, maxiter = 10)
        if i % 100 == 0:
            print i, F, 
        parameters["est_daymax"][i] = P[0]
        parameters["est_proj"][i] = P[1:]
        parameters["LC_fit_Cmat"][i] = Cmat
        
    return parameters

def M_pullfn(P, passdata):
    [parameters, data, settings] = passdata[0]

    parameters["est_EVs"] = np.reshape(P, [sum(settings["EV_treatment"] == 2), settings["nlm"], settings["nph"]])
    parameters = get_splines(parameters, settings)

    the_model = modelfn(parameters, settings, sne_to_do = None)
    pulls = np.array([], dtype=np.float64)
    for i in range(settings["nsn"]):
        these_pulls = np.sqrt(data["invvars"][i])*(the_model[i] - data["fluxes"][i])
        pulls = np.concatenate((pulls, these_pulls.flatten()))
    
    return pulls



def M_step(parameters, data, settings):
    print "Starting M step..."

    P, F, NA = miniLM_new(ministart = np.reshape(parameters["est_EVs"], sum(settings["EV_treatment"] == 2)*settings["nlm"]*settings["nph"]),
                          miniscale = np.ones(sum(settings["EV_treatment"] == 2)*settings["nlm"]*settings["nph"], dtype=np.float64),
                          residfn = M_pullfn, passdata = [parameters, data, settings], verbose = True, maxiter = 1)

    parameters["est_EVs"] = np.reshape(P, [sum(settings["EV_treatment"] == 2), settings["nlm"], settings["nph"]])
    return parameters, F

parameters, data, settings = get_settings()
parameters = get_initial_parameters(parameters, settings)

last_chi2 = 1e100
step_chi2 = 1e99

iter_count = 0

while (last_chi2 > step_chi2 + 0.1) and (iter_count < settings["max_iter"]):
    print "Starting iteration last_chi2, step_chi2", last_chi2, step_chi2, time.asctime()
    last_chi2 = step_chi2
    parameters = E_step(parameters, data, settings)
    parameters, step_chi2 = M_step(parameters, data, settings)
    parameters = get_splines(parameters, settings)
    iter_count += 1

parameters["est_EV_splines"] = None
pickle.dump([parameters, data, settings], open("results.pickle", 'wb'))
