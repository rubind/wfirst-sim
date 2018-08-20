import numpy as np
from scipy.interpolate import RectBivariateSpline
from DavidsNM import miniLM_new
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import time
import sys


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

    settings["filt_boundaries"] = 4000.*np.exp(np.arange(settings["nflt"] + 1.)*settings["filt_spacing"])
    settings["filt_lambs"] = 0.5*(settings["filt_boundaries"][1:] + settings["filt_boundaries"][:-1])
    settings["filt_names"] = [str(np.around(item))[:3] for item in settings["filt_lambs"]]
    settings["filt_colors"] = ['m', 'b', 'c', 'g',
                               'orange', 'r', 'brown', 'k']

    settings["phases"] = np.arange(-10., 36., 5.)
    settings["nph"] = len(settings["phases"])

    settings["rest_nodes"] = np.exp(np.linspace(np.log(3300.) + settings["filt_spacing"]/2.,
                                                np.log(3300.) + 1. - settings["filt_spacing"]/2.,
                                                settings["nlm"])) # Filters have log width filt_spacing

    for i in range(settings["nflt"]):
        zmin = np.clip(settings["filt_lambs"][i]/settings["rest_nodes"].max() - 1., 0., None)
        zmax = settings["filt_lambs"][i]/settings["rest_nodes"].min() - 1.

        print "%s from z=%.2f to %.2f" % (settings["filt_names"][i],
                                          zmin, zmax)

        plot_one_plus_z = np.linspace(1. + zmin, 1. + zmax, 50)

        plt.plot(plot_one_plus_z, settings["filt_lambs"][i]/plot_one_plus_z, color = settings["filt_colors"][i], linewidth = 2)
        plt.text(plot_one_plus_z[-1], settings["filt_lambs"][i]/plot_one_plus_z[-1], settings["filt_names"][i], ha = 'center', va = 'top', color = settings["filt_colors"][i])
    plt.xscale('log')

    tmp_xticks = [1., 1.3, 1.6, 2., 2.5, 3, 4., 5]
    #plt.gca().set_xticks([])
    plt.xticks(tmp_xticks, ["%.1f" % (item - 1.) for item in tmp_xticks])
    plt.savefig("filter_redshift.pdf", bbox_inches = 'tight')
    plt.close()



    print "settings:"
    for key in settings:
        print key, settings[key]


    return settings

def get_initial_parameters(settings):
    parameters = {}

    parameters["est_EVs"] = np.random.normal(size = [sum(settings["fit_ev"]), settings["nlm"], settings["nph"]])
    parameters["est_proj"] = np.ones([settings["nsn"], settings["nev"]], dtype=np.float64)
    parameters["est_daymax"] = np.zeros(settings["nsn"], dtype=np.float64)

    mean_SN_grid = np.outer(np.ones(settings["nlm"], dtype=np.float64),
                            ((settings["phases"] - - 20.)/20.)**2. * np.exp(-settings["phases"]/10.))

    parameters["est_EV_splines"] = [RectBivariateSpline(settings["rest_nodes"], settings["phases"], mean_SN_grid, kx = 3, ky = 3), None]
    parameters = get_splines(parameters, settings)

    return parameters

def get_splines(parameters, settings):
    ind = 0
    assert len(parameters["est_EVs"]) == sum(settings["fit_ev"])

    for i in range(settings["nev"]):
        if settings["fit_ev"][i]:
            parameters["est_EV_splines"][i] = RectBivariateSpline(settings["rest_nodes"], settings["phases"], parameters["est_EVs"][ind], kx = 3, ky = 3)
            ind += 1

    return parameters

def make_test_data(settings):
    true_pcs = np.random.normal(size = [settings["nsn"], 2])
    true_pcs[:,0] = 1.

    data = dict(fluxes = [], invvars = [], dates = [], redshifts = np.linspace(np.log(1.02), np.log(2.7), settings["nsn"]), true_pcs = true_pcs)


    for i in range(settings["nsn"]):
        z = data["redshifts"][i]

        dates = np.arange(-15*(1. + z), 40*(1. + z), 5.) + np.random.random()*5.
        phases = dates/(1. + z)
        dates = dates[np.where((phases >= -10) & (phases <= 35.))]
        phases = dates/(1. + z)
        
        these_fluxes = np.zeros([settings["nflt"], len(dates)], dtype=np.float64)
        these_invvars = np.zeros([settings["nflt"], len(dates)], dtype=np.float64)

        for j in range(settings["nflt"]):
            rlamb = settings["filt_lambs"][j]/(1. + z)
            if (rlamb >= settings["rest_nodes"][0]) and (rlamb <= settings["rest_nodes"][-1]):
                these_fluxes[j] = ((phases - - 20.)/20.)**2. * np.exp(-phases/10.)
                these_invvars[j] = 15**2.
                these_fluxes[j] += data["true_pcs"][i][1]*0.1*np.exp(-0.03*phases**2.)

        data["dates"].append(dates)
        data["fluxes"].append(these_fluxes)
        data["invvars"].append(these_invvars)
        
    return data

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
        rlambs = settings["filt_lambs"]/(1. + z)

        this_model = 0.
        for j in range(settings["nev"]):
            this_model += parameters["est_proj"][i][j]*parameters["est_EV_splines"][j](rlambs, phases, grid=True)

        the_model.append(this_model)
    return the_model
        

def E_pullfn(P, passdata):
    [i, parameters, data, settings] = passdata[0]

    parameters["est_daymax"][i] = P[0]
    parameters["est_proj"][i] = P[1:]

    the_model = modelfn(parameters, settings, sne_to_do = i)
    pulls = np.sqrt(data["invvars"][i])*(the_model[0] - data["fluxes"][i])
    
    return pulls.flatten()


def E_step(parameters, data, settings):
    print "Starting E step..."

    for i in range(settings["nsn"]):
        P, F, NA = miniLM_new(ministart = np.concatenate(([parameters["est_daymax"][i]], parameters["est_proj"][i])),
                              miniscale = 100*np.ones(1 + settings["nev"], dtype=np.float64),
                              residfn = E_pullfn, passdata = [i, parameters, data, settings], verbose = False, maxiter = 10)
        print i, str(P).replace("[", "").replace("]", ""), data["true_pcs"][i][1]
        parameters["est_daymax"][i] = P[0]
        parameters["est_proj"][i] = P[1:]
        
    return parameters

def M_pullfn(P, passdata):
    [parameters, data, settings] = passdata[0]

    parameters["est_EVs"] = np.reshape(P, [sum(settings["fit_ev"]), settings["nlm"], settings["nph"]])
    parameters = get_splines(parameters, settings)

    the_model = modelfn(parameters, settings, sne_to_do = None)
    pulls = np.array([], dtype=np.float64)
    for i in range(settings["nsn"]):
        these_pulls = np.sqrt(data["invvars"][i])*(the_model[i] - data["fluxes"][i])
        pulls = np.concatenate((pulls, these_pulls.flatten()))
    
    return pulls


def M_step(parameters, data, settings):
    print "Starting M step..."

    P, F, NA = miniLM_new(ministart = np.reshape(parameters["est_EVs"], sum(settings["fit_ev"])*settings["nlm"]*settings["nph"]),
                          miniscale = np.ones(sum(settings["fit_ev"])*settings["nlm"]*settings["nph"], dtype=np.float64),
                          residfn = M_pullfn, passdata = [parameters, data, settings], verbose = True, maxiter = 1)

    parameters["est_EVs"] = np.reshape(P, [sum(settings["fit_ev"]), settings["nlm"], settings["nph"]])
    return parameters, F

settings = get_settings()
parameters = get_initial_parameters(settings)
data = make_test_data(settings)

last_chi2 = 1e100
step_chi2 = 1e99

while last_chi2 > step_chi2 + 0.1:
    print "Starting iteration last_chi2, step_chi2", last_chi2, step_chi2, time.asctime()
    last_chi2 = step_chi2
    parameters = E_step(parameters, data, settings)
    parameters, step_chi2 = M_step(parameters, data, settings)
    parameters = get_splines(parameters, settings)
    
