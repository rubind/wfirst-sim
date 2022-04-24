from copy import deepcopy
from numpy import *
from DavidsNM import save_img
import time
from scipy.interpolate import interp1d
from astro_functions import CCM, get_FoM, FoM_bin_w, FoM_bin_rho, write_Cmat
import sys
#from string import strip
from FileRead import readcol
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../pixel-level/")
sys.path.append("../../pixel-level/")
from pixel_level_ETC2 import get_spec_with_err, initialize_PSFs, solve_for_exptime
plt.rcParams["font.family"] = "serif"
import os
wfirst_data_path = os.environ["WFIRST_SIM_DATA"]

####################################################### Initialize Parameters ########################################################
def get_params(the_file, PSFs = None):
    f = open(the_file)
    lines = f.read().split('\n')
    f.close()
    
    params = {"orig_lines": deepcopy(lines)}

    for i, line in enumerate(lines):
        tmp_line = line + " " # This is so that -1 goes all the way to the end of the line
        tmp_line = tmp_line[:tmp_line.find("#")]
        lines[i] = tmp_line.strip()
    
    lines = [item for item in lines if len(item) > 0]

    
    for line in lines:
        parsed = line.split(None)
        params[parsed[0]] = eval(" ".join(parsed[1:]))


    params["z_list"] = array(params["z_list"])/100.

    params["n_redshift"] = len(params["z_list"])


    min_log_restlamb = 1.e10
    max_log_restlamb = -1.e10
    
    if PSFs == None:
        PSFs = initialize_PSFs(pixel_scales = [int(around(params["pixel_scale"]/0.005))],
                               slice_scales = [int(around(params["slice_scale"]/0.005))], PSF_source = params["PSF_source"], path = wfirst_data_path + "/pixel-level/")
        
    params["weight_in_mags_perSN"] = array([], dtype=float64)
    params["redshift_vect"] = array([], dtype=float64)
    params["obs_lambs"] = array([], dtype=float64)
    params["rest_lambs"] = array([], dtype=float64)
    params["indices"] = array([], dtype=int32)

    try:
        params["exp_times"]
    except:
        params["exp_times"] = None

    params["log10_elec_per_sec"] = array([], dtype=float64)


    for i, redshift in enumerate(params["z_list"]):
        if i == 0:
            print("Working on nearby SNe")
            
            obs_waves = exp(arange(log(params["min_ground_wavelength"]), log(params["max_ground_wavelength"]), 0.01))

            good_inds = where((obs_waves/(1. + redshift) >= params["min_rest_wavelength"])*(obs_waves/(1. + redshift) <= params["max_rest_wavelength"]))
            nwave = len(obs_waves[good_inds])
            params["weight_in_mags_perSN"] = concatenate((params["weight_in_mags_perSN"], ones(nwave, dtype=float64)*100.))
            params["log10_elec_per_sec"] = concatenate((params["log10_elec_per_sec"], zeros(nwave, dtype=float64)))
            params["redshift_vect"] = concatenate((params["redshift_vect"], ones(nwave, dtype=float64)*redshift))
            params["obs_lambs"] = concatenate((params["obs_lambs"], obs_waves[good_inds]))
            params["rest_lambs"] = concatenate((params["rest_lambs"], obs_waves[good_inds]/(1. + redshift)))
            params["indices"] = concatenate((params["indices"], zeros(nwave, dtype=int32)))

            nearby_SN_e_per_sec = get_spec_with_err(redshift, exp_time = 1., phase = 0, gal_flamb = lambda x:0., pixel_scale = params["pixel_scale"], slice_scale = params["slice_scale"],
                                                    show_plots = 0, IFURfl = params["IFURfl"], min_wave = params["min_wavelength"], max_wave = params["max_wavelength"], PSFs = PSFs,
                                                    source_dir = wfirst_data_path + "/pixel-level/input/")["PSF_wghtd_e_SN"]

        else:

            if params["exp_times"] == None or params["exp_times"] == []:
                exp_time = solve_for_exptime(10., redshift, PSFs, key1 = "rest_frame_band_S/N", key2 = (5000, 6000),
                                             pixel_scale = params["pixel_scale"], slice_scale = params["slice_scale"],
                                             source_dir = wfirst_data_path + "/pixel-level/input/",
                                             IFURfl = params["IFURfl"], min_wave = params["min_wavelength"], max_wave = params["max_wavelength"])
            else:
                exp_time = params["exp_times"][i-1]
            
            signal_to_noises = get_spec_with_err(redshift, exp_time, phase = 0, gal_flamb = lambda x:0., pixel_scale = params["pixel_scale"], slice_scale = params["slice_scale"],
                                                 show_plots = 0, IFURfl = params["IFURfl"], min_wave = params["min_wavelength"], max_wave = params["max_wavelength"], PSFs = PSFs,
                                                 source_dir = wfirst_data_path + "/pixel-level/input/")
            
            
            good_inds = where((signal_to_noises["obs_waves"]/(1. + redshift) >= params["min_rest_wavelength"])*(signal_to_noises["obs_waves"]/(1. + redshift) <= params["max_rest_wavelength"]))
            nwave = len(good_inds[0])


            params["log10_elec_per_sec"] = concatenate((params["log10_elec_per_sec"], log10(signal_to_noises["PSF_wghtd_e_SN"]/(exp_time*nearby_SN_e_per_sec))[good_inds]
                                                    ))

            params["weight_in_mags_perSN"] = concatenate((params["weight_in_mags_perSN"], signal_to_noises["spec_S/N"][good_inds]**2. * 0.848304))
            params["redshift_vect"] = concatenate((params["redshift_vect"], ones(nwave, dtype=float64)*redshift))
            params["obs_lambs"] = concatenate((params["obs_lambs"], signal_to_noises["obs_waves"][good_inds]))
            params["rest_lambs"] = concatenate((params["rest_lambs"], signal_to_noises["obs_waves"][good_inds]/(1. + redshift)))
            params["indices"] = concatenate((params["indices"], ones(nwave, dtype=int32)*i))


    #plt.plot(params["log10_elec_per_sec"])
    #plt.show()

    #plt.plot(params["obs_lambs"], params["weight_in_mags"], '.')
    #plt.show()

    params["model_rest_lambs"] = exp(arange(log(min(params["rest_lambs"])) - 0.006, log(max(params["rest_lambs"])) + 0.01, 0.01))
    params["model_obs_lambs"] = exp(arange(log(min(params["obs_lambs"])) - 0.006, log(max(params["obs_lambs"])) + 0.01, 0.01))

    params["n_obslambs"] = len(params["model_obs_lambs"])
    params["n_restlambs"] = len(params["model_rest_lambs"])
    params["n_dex"] = int(around(ceil(max(   abs(params["log10_elec_per_sec"])   ))))

    print("params read ", params)

    return params

####################################################### Parse Vector to Dict ########################################################


def parseP(P, params):
    parsed = {}
    ind = 0
    
    nred = params["n_redshift"]
    nlamb = params["n_obslambs"]
    nrestlamb = params["n_restlambs"]

    parsed["mag"] = P[ind:ind+nred]
    ind += nred

    parsed["mag_nuisance"] = P[ind:ind+nred]
    ind += nred

    parsed["EBV"] = P[ind:ind+nred]
    ind += nred

    for i, spectral_feature_name in enumerate(params["spectral_feature_names"]):
        parsed[spectral_feature_name] = P[ind:ind+3]
        ind += 3

    parsed["mean"] = P[ind:ind+nrestlamb]
    ind += nrestlamb

    parsed["ZP_indep"] = P[ind:ind+nlamb]
    ind += nlamb
    parsed["ZP_slope"] = P[ind]
    ind += 1
    parsed["ZP_ground_to_space"] = P[ind]
    ind += 1

    parsed["ZP_dex"] = P[ind:ind + params["n_dex"] + 1]
    ind += params["n_dex"] + 1

    return parsed

####################################################### Parse Dict to Vector  ########################################################

def unparseP(parsed, params):
    P = concatenate((parsed["mag"], parsed["mag_nuisance"], parsed["EBV"]))
    
    for spectral_feature_name in params["spectral_feature_names"]:
        P = append(P, parsed[spectral_feature_name])

    P = concatenate((P, parsed["mean"], parsed["ZP_indep"], [parsed["ZP_slope"], parsed["ZP_ground_to_space"]], parsed["ZP_dex"]))
    return P


####################################################### Make the Model  ########################################################

def get_bin_value(the_array, redshifts, i):
    if len(the_array) == len(redshifts):
        # One value for each bin
        if i == None:
            return the_array
        else:
            return the_array[i]
    else:
        # Quadratic in redshift, specified at z=0, 0.5, 1.0
        quad_a = 2*the_array[0] - 4.*the_array[1] + 2*the_array[2]
        quad_b = -3*the_array[0] + 4.*the_array[1] - the_array[2]
        quad_c = the_array[0]

        if i == None:
            return quad_c + quad_b*redshifts + quad_a*redshifts**2.
        else:
            return quad_c + quad_b*redshifts[i] + quad_a*redshifts[i]**2.


def modelfn(parsed, params):
    the_models = array([], dtype=float64) # Each redshift bin, in turn
    HD_RMSs = array([], dtype=float64)

    for i, redshift in enumerate(params["z_list"]):
        HD_RMSs = append(HD_RMSs, sqrt(params["HD_RMS"][i]**2. + (params["HD_lensing"]*redshift)**2. + (5*params["v_pec"]/(log(10.)*299792.458*redshift))**2.
                                   ))
    

    model = parsed["mag"][params["indices"]] + parsed["mag_nuisance"][params["indices"]]
    model += CCM(params["rest_lambs"], R_V = params["R_V"])*parsed["EBV"][params["indices"]]*params["R_V"]
    model += interp1d(params["model_rest_lambs"], parsed["mean"], kind = 'linear')(params["rest_lambs"])

    model += interp1d(params["model_obs_lambs"], parsed["ZP_indep"], kind = 'linear')(params["obs_lambs"])
    model += interp1d(arange(0., -params["n_dex"] - 0.01, -1.), parsed["ZP_dex"], kind = 'linear')(params["log10_elec_per_sec"])
    model += parsed["ZP_slope"]*(params["obs_lambs"] - 10000.)/10000. # Slope per micron
        
    lowzinds = where(params["indices"] == 0)
    model[lowzinds] += 0.5*parsed["ZP_ground_to_space"]
    highzinds = where(params["indices"] != 0)
    model[highzinds] -= 0.5*parsed["ZP_ground_to_space"]


    for i in range(len(params["z_list"])):
        for this_name, this_wavelength, this_width, this_slope in zip(params["spectral_feature_names"], params["spectral_feature_wavelengths"], params["spectral_feature_widths"], params["spectral_feature_correction_slopes"]):
            this_depth = get_bin_value(parsed[this_name], params["z_list"], i)
        
            inds = where(params["indices"] == i)
            model[inds] += this_depth*(exp(-0.5*(params["rest_lambs"][inds] - this_wavelength)**2. / this_width**2.) + this_slope)

    return model, HD_RMSs


def residfn(P, wrapped):
    params = wrapped[0]

    parsed = parseP(P, params)
    model, HD_RMSs = modelfn(parsed, params)
    #print model


    pulls = model*sqrt(params["weight_in_mags_perSN"])
    assert all(1 - isnan(pulls)), "1"

    pulls = append(pulls, parsed["mag_nuisance"]/HD_RMSs)
    assert all(1 - isnan(pulls)), "2"

    pulls = append(pulls, parsed["ZP_indep"]/params["ZP_indep"])
    assert all(1 - isnan(pulls)), "3"

    pulls = append(pulls, parsed["ZP_slope"]/params["ZP_slope"])
    assert all(1 - isnan(pulls)), "4"


    pulls = append(pulls, parsed["ZP_ground_to_space"]/params["ZP_ground_to_space"])

    assert all(1 - isnan(pulls)), "5"

    pulls = append(pulls, parsed["ZP_dex"]/params["ZP_dex"])

    pulls = append(pulls, parsed["mean"].sum() / params["mean_norm"])
    pulls = append(pulls, dot(parsed["mean"], CCM(params["model_rest_lambs"], R_V = params["R_V"])) / params["EBV_norm"])

    
    for this_wavelength, this_width in zip(params["spectral_feature_wavelengths"], params["spectral_feature_widths"]):
        pulls = append(pulls, dot(exp(-0.5*(params["model_rest_lambs"] - this_wavelength)**2. / this_width**2.), parsed["mean"])/params["feature_norm"])

    #print "pull"
    return pulls
"""
def plot_summ():
    plt.figure(figsize = (8, len(parsed_errs.keys())*3))

    for i, key in enumerate(parsed_errs):
        plt.subplot(len(parsed_errs.keys()), 1, i + 1)


        y_scale = 1
        if type(parsed_errs[key]) == type(arange(2.)[0]):
            x_vals = 0
        elif len(parsed_errs[key]) == params["n_redshift"]:
            x_vals = params["z_list"]
            y_scale = sqrt(array(params["NSNe"]))
        elif len(parsed_errs[key]) == params["n_restlambs"]:
            x_vals = params["model_rest_lambs"]
        elif len(parsed_errs[key]) == params["n_obslambs"]:
            x_vals = params["model_obs_lambs"]
        elif len(parsed_errs[key]) == 3:
            x_vals = [0, 0.5, 1.]
        else:
            print key
            sys.exit(1)

        plt.plot(x_vals, parsed_errs[key]*y_scale, 'o')
        if key == "mag":
            plt.plot(x_vals, parsed_errs["mag_nuisance"]*y_scale, 'o')
            plt.plot(x_vals, y_scale*sqrt(parsed_errs[key]**2. - parsed_errs["mag_nuisance"]**2.))
        ylim = plt.ylim()
        plt.ylim(0, ylim[1])
        plt.title(key, size = 7)
        plt.xticks(size = 7)
        plt.yticks(size = 7)
    plt.savefig("errs.pdf")

    cor_mat = zeros(Cmat.shape, dtype=float64)
    for i in range(len(Cmat)):
        for j in range(i, len(Cmat)):
            cor_mat[i,j] = Cmat[i,j]/sqrt(Cmat[i,i]*Cmat[j,j])
            cor_mat[j,i] = cor_mat[i,j]

    save_img(cor_mat, "cor_mat.fits")

    #plt.figure(figsize = (8,4))
    #plt.subplot(1,2,1)
    plt.imshow(sn_Cmat, interpolation = 'nearest', origin = 'lower', cmap=plt.get_cmap('Greys'), aspect = 'equal')
    #plt.colorbar()
    plt.xticks([0,5,10,15], ["0.05", "0.55", "1.05", "1.55"])
    #plt.xticklabels(["0.05", "0.55", "1.05", "1.55"])
    plt.yticks([0,5,10,15], ["0.05", "0.55", "1.05", "1.55"])
    #plt.yicklabels(["0.05", "0.55", "1.05", "1.55"])
    #plt.subplot(1,2,2)
    #plt.plot(params["z_list"], sqrt(diag(sn_Cmat)*params["NSNe"]), color = 'k')
    #plt.xlabel("Redshift")
    #plt.ylabel("Total Distance Modulus Uncertainty per SN")
    plt.title("Estimated Covariance Matrix, Binned in Redshift")
    plt.xlabel("Redshift")
    plt.ylabel("Redshift")

    plt.savefig("Cmat.eps", bbox_inches = 'tight')

    write_Cmat(params["orig_lines"], sn_Cmat, FoM, flname = "Cmat_" + params["suffix"] + ".txt")
"""

def get_Jacobian_NSNe1(params):
    ministart = {"mag": [0]*params["n_redshift"], "mag_nuisance": [0]*params["n_redshift"],
                 "EBV": [0]*params["n_redshift"], "mean": [0]*params["n_restlambs"],
                 "ZP_indep": [0]*params["n_obslambs"], "ZP_slope": 0, "ZP_ground_to_space": 0,
                 "ZP_dex": [0]*(params["n_dex"] + 1)}

    
    for spectral_feature_name in params["spectral_feature_names"]:
        ministart[spectral_feature_name] = [0,0,0]

    print(ministart)
        


    ministart = array(unparseP(ministart, params), dtype=float64)
    print(ministart)


    orig = residfn(ministart, [params])
    miniscale = 1e-6

    jacobian = zeros([len(orig), len(ministart)], dtype=float64)
    print(jacobian.shape)
    
    for i in range(len(ministart)):
        new_pos = deepcopy(ministart)
        new_pos[i] += miniscale

        new = residfn(new_pos, [params])
        jacobian[:,i] = (new - orig)/miniscale
    return jacobian

def get_Jacobian(jacobian_NSNe1, params):
    jacobian = deepcopy(jacobian_NSNe1)

    NSNe_by_wave = zeros(len(params["indices"]), dtype=float64)
    
    for i in range(len(params["NSNe"])):
        NSNe_by_wave += params["NSNe"][i]*(params["indices"] == i)
        
    for i in range(len(NSNe_by_wave)):
        jacobian[i] *= sqrt(NSNe_by_wave[i])

    for i in range(len(params["NSNe"])):
        jacobian[len(NSNe_by_wave) + i] *= sqrt(params["NSNe"][i])
    return jacobian



def run_FoM(paramfl = None, PSFs = None, jacobian_NSNe1 = None, params = None):
    if params == None:
        params = get_params(paramfl, PSFs)

    show_plots = 0

    if all(jacobian_NSNe1 == None):
        jacobian_NSNe1 = get_Jacobian_NSNe1(params)

    jacobian = get_Jacobian(jacobian_NSNe1, params)

    print(time.time())
    Wmat = dot(transpose(jacobian), jacobian)
    print(time.time())
    Cmat = linalg.inv(Wmat)
    print(time.time())

    if show_plots:
        save_img(Cmat, "Cmat.fits")

    errs = sqrt(diag(Cmat))

    parsed_errs = parseP(errs, params)


    test_inv = dot(Cmat, linalg.inv(Cmat))
    print("Inverse Test", abs(test_inv - identity(len(test_inv))).max())
    assert abs(test_inv - identity(len(test_inv))).max() < 1.e-6, "Inversion problem!"

    sn_Cmat = Cmat[:len(params["z_list"]), :len(params["z_list"])]
    #save_img(sn_Cmat, "sn_Cmat.fits")

    if params["FoM_type"][0] == "DETF":
        FoM, uncertainties = get_FoM(sn_Cmat, params["z_list"], zp = 0.3)
    elif params["FoM_type"][0] == "Binw":
        FoM, uncertainties = FoM_bin_w(params["z_list"], sn_Cmat, bins = params["FoM_type"][1:])
    elif params["FoM_type"][0] == "Binrho":
        FoM, uncertainties = FoM_bin_rho(params["z_list"], sn_Cmat, bins = params["FoM_type"][1:])
    print("FoM", FoM)
    return FoM


####################################################### Main ########################################################

if __name__ == "__main__":
    run_FoM(sys.argv[1])
