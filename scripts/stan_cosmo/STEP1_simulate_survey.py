from numpy import *
from astropy.cosmology import FlatLambdaCDM
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.interpolate import interp1d
import sys
import os
wfirst_path = os.environ["WFIRST"]
wfirst_data_path = os.environ["WFIRST_SIM_DATA"]
sys.path.append(wfirst_path + "/scripts/pixel-level/")
sys.path.append(wfirst_path + "/scripts/host/")
from host_generator import make_galaxy_spectrum
from generator import make_SALT2_params
from pixel_level_ETC2 import initialize_PSFs, get_imaging_SN, solve_for_exptime, get_spec_with_err
#from FileRead import readcol
import subprocess
import argparse
import glob
import sncosmo
import pickle as pickle
import multiprocessing as mp
from astropy.table import Table, vstack
import matplotlib.path as MPLpath
import copy as cp
import tempfile
import time
from scipy.stats import percentileofscore
import matplotlib.patches as patches
whoami = subprocess.getoutput("whoami")


def file_to_fn(fl, col = 1):
    vals = loadtxt(fl)
    x = vals[:,0]
    y = vals[:,col]

    return interp1d(x, y, kind = 'linear', bounds_error = False, fill_value = 0.)

"""
def eval_file(param):
    #f = open(opts.fl)
    f = open(sys.argv[1])
    lines = f.read().split('\n')
    f.close()

    for line in lines:
        parsed = line.split(None)
        if len(parsed) > 1:
            if parsed[0] == param:
                return eval(" ".join(parsed[1:]))
    print "Couldn't find ", param
    raise Exception("Missing Value")
"""

def read_from_lines(lines, param, islist = False, item_number = 1):
    found_param = 0

    for line in lines:
        parsed = line.split(',')
        if len(parsed) > 1:
            if parsed[0] == param:
                found_param += 1

                if found_param == item_number:
                    vals_to_return = []
                    for item in parsed[1:]:
                        if item != "":
                            try:
                                vals_to_return.append(eval(item))
                            except:
                                vals_to_return.append(item)
                    if islist:
                        return vals_to_return
                    else:
                        return vals_to_return[0]
    return None

def read_csv(csv_file):
    f = open(csv_file)
    lines = f.read().replace('\r', '\n').split('\n')
    f.close()

    survey_parameters = {}

    # Start with global parameters
    for key in ["total_survey_time", "survey_duration", "hours_per_visit", "maximum_trigger_fraction", "shallow_SNR", "medium_SNR", "deep_SNR", "reference_SNR", "number_of_reference_dithers", "SN_rates",
                "slew_table", "zodiacal_background", "telescope_temperature", "PSFs", "WFI_PSFs", "interpixel_capacitance",
                "WFI_dark_current", "WFI_read_noise_floor", "WFI_read_noise_white", "WFI_pixel_scale",
                "IFU_min_wave", "IFU_max_wave", "IFU_effective_area", "IFU_resolution", "IFU_pixel_scale", "IFU_slice_in_pixels", "IFU_dark_current", "IFU_read_noise_floor", "IFU_read_noise_white", "bad_pixel_rate"]:
        survey_parameters[key] = read_from_lines(lines, key)
    for key in ["adjust_each_SN_exp_time", "targeted_parallels"]:
        survey_parameters[key] = ["true", "t", "yes", "y"].count(read_from_lines(lines, key).lower())
    for key in ["normalization_wavelength_range", "spectra_depth_and_phase", "grizY_30s_ground_depths"]:
        survey_parameters[key] = read_from_lines(lines, key, islist = True)
    survey_parameters["spectra_depth_and_phase"] = [(item.split(None)[0], float(item.split(None)[1])) for item in survey_parameters["spectra_depth_and_phase"]]

    # Now read each tier
    n_tiers = 0
    while read_from_lines(lines, "tier_name", item_number = n_tiers + 1) != None:
        n_tiers += 1
    print("n_tiers ", n_tiers)
    

    single_keys = ["tier_name", "tier_fraction_time", "square_degrees", "cadence"]
    list_keys = ["filters", "exp_times_per_dither", "dithers_per_filter", "trigger_redshifts", "trigger_fraction", "parallel_filters"]
    survey_parameters["tier_parameters"] = {}
    for key in single_keys + list_keys:
        survey_parameters["tier_parameters"][key] = []


    for i in range(n_tiers):
        for key in single_keys:
            survey_parameters["tier_parameters"][key].append(read_from_lines(lines, key, item_number = i + 1, islist = False))
        for key in list_keys:
            survey_parameters["tier_parameters"][key].append(read_from_lines(lines, key, item_number = i + 1, islist = True))

    print(survey_parameters)
    return survey_parameters



def get_ground_AB_mag(flamb_SN, filt):
    AB_mag = -2.5*log10(sum(ground_filt_fns[filt](ground_obslambs)*flamb_SN*ground_obslambs)/
                        sum(ground_filt_fns[filt](ground_obslambs)*(0.108848062485/ground_obslambs**2.)*ground_obslambs))
    return AB_mag

def init_ground(grizY_30s_ground_depths):
    ground_obslambs = arange(3000., 11000., 10.)
    ground_filt_fns = {}

    for filt in "grizY":
        ground_filt_fns[filt] = file_to_fn(wfirst_data_path + "/pixel-level/input/LSST_" + filt + ".txt")

    ground_five_sigma_one_hour = dict(g = grizY_30s_ground_depths[0] + 1.25*log10(3600/30.),
                                      r = grizY_30s_ground_depths[1] + 1.25*log10(3600/30.),
                                      i = grizY_30s_ground_depths[2] + 1.25*log10(3600/30.),
                                      z = grizY_30s_ground_depths[3] + 1.25*log10(3600/30.),
                                      Y = grizY_30s_ground_depths[4] + 1.25*log10(3600/30.))

    WFI_filt_fns = {}
    for filt in ["R062", "Z087", "Y106", "J129", "H158", "F184", "K213", "W146", "Euclid_Y", "Euclid_J", "Euclid_H"]:
        WFI_filt_fns[filt] = file_to_fn(wfirst_data_path + "/pixel-level/input/" + filt + ".txt")

    return ground_filt_fns, ground_obslambs, ground_five_sigma_one_hour, WFI_filt_fns

def get_ground_depths(survey_parameters, tier):
    these_filts = survey_parameters["tier_parameters"]["filters"][tier]
    these_exps = survey_parameters["tier_parameters"]["exp_times_per_dither"][tier]

    ground_depths = {}
    for j in range(len(these_filts)):
        if len(these_filts[j]) == 1:
            # Ground filter found
            print("Ground filter found ", these_filts[j])
            assert survey_parameters["tier_parameters"]["dithers_per_filter"][tier][j] == 1, "It's a waste to model more than one ground dither!"
            ground_depths[these_filts[j]] = ground_five_sigma_one_hour[these_filts[j]] + 1.25*log10(these_exps[j]/3600.)
    print("ground_depths found ", ground_depths)
    return ground_depths
    

def get_SNCosmo_model(redshift, x1, c, MV, daymax, source):
    #sncosmo_model = sncosmo.Model(source="salt2-extended")
    sncosmo_model = sncosmo.Model(source=source)

    cosmo = FlatLambdaCDM(Om0 = 0.3, H0 = 70.)
    ten_pc_z = 2.33494867e-9
    assert abs(cosmo.distmod(z=ten_pc_z).value) < 1.e-3, "Distance modulus zeropoint wrong!"

    ampl = 1.
    for i in range(2):
        sncosmo_model.set(z=ten_pc_z, t0=0., x0=ampl, x1 = x1, c = c)

        mag = sncosmo_model.bandmag('bessellv', 'ab', 0.)
        #print "mag, ampl ", mag, ampl
        ampl *= 10.**(0.4*(mag - MV))
        
    mu = cosmo.distmod(redshift).value
    #print "mu ", mu

    sncosmo_model.set(z=redshift, t0=daymax, x0 = ampl*10**(-0.4*mu))
    return sncosmo_model


def realize_SN(redshift, daymax, source):
    """
    x1 = random.normal() - 0.25
    c = random.normal()*0.1 + random.exponential()*0.1 - 0.1*log(2.)
    host_mass = random.normal() + 10.

    MV = -19.08 + random.normal()*0.12 + random.normal()*0.055*redshift - 0.13*x1 + 2.1*c - (host_mass > 10.)*0.08 # beta is 2.1, as this is MV
    """
    
    MB, x1, color, mass = make_SALT2_params(size = 1)
    [MB, x1, color, mass] = [item[0] for item in (MB, x1, color, mass)]

    MB += random.normal()*0.1 + 0.055*redshift*random.normal() + 5/log(10.)*(0.001/redshift)*random.normal()

    sncosmo_model = get_SNCosmo_model(redshift = redshift,
                                      x1 = x1,
                                      c = color,
                                      MV = MB - color,
                                      daymax = daymax, source = source)


    return MB - color, x1, color, mass, sncosmo_model


def round_to_cadence(vals, cadence):
    return cadence*around(vals/float(cadence))

def date_to_roll_angle(the_date):
    even_step = float(isclose(the_date % 2, 0)) # Assumes cadence is odd!!!!!!
    delta_angle = even_step*0.35 - 0.175
    return (the_date*2.*pi/365.24 + delta_angle) % (2.*pi)

def compress_lc_data(SN_data, current_date, filt):
    """If more than one observation per day, bin."""

    

    for i in range(SN_data["nsne"]):
        inds = where(isclose(SN_data["SN_observations"][i]["dates"], current_date, rtol = 0)*
                     (SN_data["SN_observations"][i]["filts"] == filt))[0]
        
        if len(inds) > 1:
            # Time to compress
            
            nobs = len(SN_data["SN_observations"][i]["dfluxes"])

            assert all(inds == arange(nobs - len(inds), nobs)) # These should be the last observations added
            
            ws = SN_data["SN_observations"][i]["dfluxes"][inds]**(-2.)
            fs = SN_data["SN_observations"][i]["fluxes"][inds]

            new_f = sum(ws * fs)/sum(ws)
            new_df = 1./sqrt(sum(ws))

            SN_data["SN_observations"][i]["fluxes"][inds[0]] = new_f
            SN_data["SN_observations"][i]["dfluxes"][inds[0]] = new_df
            SN_data["SN_observations"][i]["ispar"][inds[0]] = any(SN_data["SN_observations"][i]["ispar"][inds])

            for key in ["dates", "fluxes", "dfluxes", "filts"]:
                SN_data["SN_observations"][i][key] = SN_data["SN_observations"][i][key][:nobs+1-len(inds)]

    return SN_data

"""
SN_data = {"nsne": 1, "SN_observations":
           [{"dates": array([1, 5, 5, 5]),
             "filts": array(["H158", "H158", "J129", "J129"]),
             "fluxes": arange(4.), "dfluxes": array([1., 1., 1., 2.])}]}

print compress_lc_data(SN_data, 5., "J129")
"""

def field_outline_old(WFI_RA, WFI_dec, WFI_orient):
    WFI_dRA = sqrt(0.28*0.5)*0.5 # x
    WFI_ddec = sqrt(0.28*2.)*0.5 # y

    assert WFI_orient <= 2.*pi + 1e-6, "angles should be in radians!"

    x1 = WFI_RA + WFI_dRA*cos(WFI_orient) - WFI_ddec*sin(WFI_orient)
    y1 = WFI_dec + WFI_ddec*cos(WFI_orient) + WFI_dRA*sin(WFI_orient) 

    x2 = WFI_RA - WFI_dRA*cos(WFI_orient) - WFI_ddec*sin(WFI_orient)
    y2 = WFI_dec + WFI_ddec*cos(WFI_orient) - WFI_dRA*sin(WFI_orient) 

    x3 = WFI_RA - WFI_dRA*cos(WFI_orient) + WFI_ddec*sin(WFI_orient)
    y3 = WFI_dec - WFI_ddec*cos(WFI_orient) - WFI_dRA*sin(WFI_orient) 

    x4 = WFI_RA + WFI_dRA*cos(WFI_orient) + WFI_ddec*sin(WFI_orient)
    y4 = WFI_dec - WFI_ddec*cos(WFI_orient) + WFI_dRA*sin(WFI_orient) 
    return (x1, x2, x3, x4), (y1, y2, y3, y4)

def inside_field_old(SN_RAs, SN_decs, WFI_RA, WFI_dec, WFI_orient):
    (x1, x2, x3, x4), (y1, y2, y3, y4) = field_outline(WFI_RA, WFI_dec, WFI_orient)

    #plt.plot(WFI_RA, WFI_dec, 'o')
    #plt.plot([x1,x2,x3,x4,x1], [y1, y2, y3, y4, y1])

    outline = MPLpath.Path(array([[x1, y1],
                                  [x2, y2],
                                  [x3, y3],
                                  [x4, y4]]
                             ))
    try:
        return outline.contains_points(list(zip(SN_RAs, SN_decs)))
    except:
        print("Couldn't run contains_points ", SN_RAs, SN_decs)
        sys.exit(1)


def field_outline(WFI_RA, WFI_dec, WFI_orient):
    assert WFI_orient <= 2.*pi + 1e-6, "angles should be in radians!"

    WFI_dRAs = array([[-0.21786, -0.09428, -0.09461, -0.21812], [-0.08827,  0.03332,  0.03291, -0.08862], [ 0.03921,  0.15842,  0.15794,  0.0388 ], [-0.25473, -0.13077, -0.13175, -0.25547], [-0.12475, -0.00264, -0.00387, -0.12574], [ 0.00328,  0.12311,  0.12167,  0.00204], [-0.3105 , -0.18618, -0.18764, -0.31155], [-0.18013, -0.05747, -0.05934, -0.18161], [-0.05152,  0.06906,  0.06681, -0.0534 ], [-0.09461, -0.21812, -0.21786, -0.09428], [ 0.03291, -0.08862, -0.08827,  0.03332], [ 0.15794,  0.0388 ,  0.03921,  0.15842], [-0.13175, -0.25547, -0.25473, -0.13077], [-0.00387, -0.12574, -0.12475, -0.00265], [ 0.12166,  0.00204,  0.00328,  0.1231 ], [-0.18764, -0.31155, -0.3105 , -0.18618], [-0.05934, -0.18161, -0.18013, -0.05747], [ 0.06681, -0.0534 , -0.05152,  0.06906]])

    WFI_ddecs = array([[-0.00308, -0.00307, -0.12839, -0.12898], [-0.00307, -0.00305, -0.12762, -0.12836], [-0.00305, -0.00303, -0.1267, -0.12758], [-0.14452, -0.14391, -0.26885, -0.26998], [-0.14387, -0.14309, -0.26733, -0.26879], [-0.14305, -0.14211, -0.26551, -0.26725], [-0.28569, -0.28464, -0.40868, -0.41017], [-0.28458, -0.28318, -0.40662, -0.4086], [-0.28311, -0.28139, -0.40407, -0.40651], [0.12839, 0.12898, 0.00308, 0.00307], [0.12762, 0.12836, 0.00307, 0.00305], [0.1267, 0.12758, 0.00305, 0.00303], [0.26885, 0.26998, 0.14452, 0.14391], [0.26733, 0.26878, 0.14387, 0.14309], [0.2655, 0.26725, 0.14305, 0.14211], [0.40868, 0.41017, 0.28568, 0.28463], [0.40661, 0.40859, 0.28457, 0.28318], [0.40406, 0.4065, 0.28311, 0.28139]])

    WFI_dRAs -= mean(WFI_dRAs)
    WFI_ddecs -= mean(WFI_ddecs)

    xs = WFI_RA + WFI_dRAs*cos(WFI_orient) - WFI_ddecs*sin(WFI_orient)
    ys = WFI_dec + WFI_ddecs*cos(WFI_orient) + WFI_dRAs*sin(WFI_orient)

    return xs, ys

"""
xs, ys = field_outline(0, 0, 0.1)
xs = concatenate((xs, xs[:,:1]), axis = 1)
ys = concatenate((ys, ys[:,:1]), axis = 1)

for i in range(len(xs)):
    plt.plot(xs[i], ys[i])


xs, ys = field_outline(0, 0, 0.3)
xs = concatenate((xs, xs[:,:1]), axis = 1)
ys = concatenate((ys, ys[:,:1]), axis = 1)

for i in range(len(xs)):
    plt.plot(xs[i], ys[i])


plt.axes().set_aspect('equal')
plt.savefig("test.pdf")
plt.close()
stop_here
"""


def inside_field(SN_RAs, SN_decs, WFI_RA, WFI_dec, WFI_orient):
    xs, ys = field_outline(WFI_RA, WFI_dec, WFI_orient)

    #plt.plot(WFI_RA, WFI_dec, 'o')
    #plt.plot([x1,x2,x3,x4,x1], [y1, y2, y3, y4, y1])

    inside_mask = 0
    
    for i in range(len(xs)):
        outline = MPLpath.Path(list(zip(xs[i], ys[i])))
        try:
            inside_mask += outline.contains_points(list(zip(SN_RAs, SN_decs)))
        except:
            print("Couldn't run contains_points ", SN_RAs, SN_decs)
            sys.exit(1)
    return inside_mask


"""
SN_RAs = random.normal(size = 10000)*0.25
SN_decs = random.normal(size = 10000)*0.25

mask = inside_field(SN_RAs = SN_RAs, SN_decs = SN_decs, WFI_RA = 0.0, WFI_dec = 0.0, WFI_orient = 0)
plt.plot(SN_RAs[where(mask)], SN_decs[where(mask)], '.', markersize = 1)
plt.plot(SN_RAs[where(1 - mask)], SN_decs[where(1 - mask)], '.', color = 'r', markersize = 1)

plt.axes().set_aspect('equal', 'datalim')
plt.savefig("inside_test.pdf")
plt.close()
stop_here
"""

def IFS_to_WFI(IFS_RA, IFS_Dec, WFI_orient):
    dr = 0.491 + 0.0712295833333 + 0.039
    return IFS_RA + cos(WFI_orient)*dr, IFS_Dec + sin(WFI_orient)*dr



def run_observation_through_ETC(SN_data, row_to_add, current_date):
    phases = (current_date - SN_data["SN_table"]["daymaxes"])/(1. + SN_data["SN_table"]["redshifts"])

    active_SN_mask = (phases >= -15.)*(phases <= 45.)
    if any(active_SN_mask):
        inds = where(active_SN_mask)
        
        #print "inds ", inds

        in_pointing_mask = inside_field(SN_RAs = SN_data["SN_table"]["RAs"][inds],
                                        SN_decs = SN_data["SN_table"]["Decs"][inds],
                                        WFI_RA = row_to_add["RA"],
                                        WFI_dec = row_to_add["dec"],
                                        WFI_orient = row_to_add["orient"])
        
        #print "in_pointing_mask ", in_pointing_mask
        #print "where(in_pointing_mask)", where(in_pointing_mask)
        inds = array(inds[0])[where(in_pointing_mask)]
        #print "inds ", inds

        for ind in inds:
            f_lamb_SN = SN_data["SN_observations"][ind]["sncosmo_model"].flux(current_date, WFI_args["waves"])

            #args = cp.deepcopy(WFI_args)
            WFI_args["mdl"] = f_lamb_SN

            ETC_result = get_imaging_SN(redshift = 0, exp_time = row_to_add["exptime"], gal_flamb = SN_data["SN_observations"][ind]["gal_background"],
                                        effective_meters2_fl = WFI_filt_fns[row_to_add["filt"]], phase = 0,
                                        offset_par = random.random()*22, offset_perp = random.random()*22, **WFI_args)

            #ETC_result = {"AB_mag": 24, "PSF_phot_S/N": 5.}

            AB_flux = 10.**(-0.4*(ETC_result["AB_mag"] - master_zp))
            SNR_total = ETC_result["PSF_phot_S/N"]

            flux_with_noise = AB_flux + random.normal()*AB_flux/SNR_total
            SN_data["SN_observations"][ind]["dates"] = append(SN_data["SN_observations"][ind]["dates"], current_date)
            SN_data["SN_observations"][ind]["fluxes"] = append(SN_data["SN_observations"][ind]["fluxes"], flux_with_noise)
            SN_data["SN_observations"][ind]["dfluxes"] = append(SN_data["SN_observations"][ind]["dfluxes"], AB_flux/SNR_total)
            SN_data["SN_observations"][ind]["filts"] = append(SN_data["SN_observations"][ind]["filts"], row_to_add["filt"])
            SN_data["SN_observations"][ind]["ispar"] = append(SN_data["SN_observations"][ind]["ispar"], row_to_add["SNind"] > -1)


    in_pointing_mask = inside_field(SN_RAs = SN_data["test_points"][0]["RAs"],
                                    SN_decs = SN_data["test_points"][0]["Decs"],
                                    WFI_RA = row_to_add["RA"],
                                    WFI_dec = row_to_add["dec"],
                                    WFI_orient = row_to_add["orient"])
    
    inds = where(in_pointing_mask)
    for ind in inds[0]:
        f_lamb_SN = 0.10884806248*10.**(-0.4*SN_data["test_points"][0]["Mags"][ind]) /(WFI_args["waves"]**2.)

        WFI_args["mdl"] = f_lamb_SN

        ETC_result = get_imaging_SN(redshift = 0, exp_time = row_to_add["exptime"], gal_flamb = lambda x:0,
                                    effective_meters2_fl = WFI_filt_fns[row_to_add["filt"]], phase = 0,
                                    offset_par = random.random()*22, offset_perp = random.random()*22, **WFI_args)

        SNR_total = ETC_result["PSF_phot_S/N"]
        SN_data["test_points"][0]["Observations"][ind][row_to_add["filt"]].append(SNR_total)
    

    if row_to_add["SNind"] > -1:
        # If IFS observation:
        ind = row_to_add["SNind"]

        if phases[ind] < 40:
            f_lamb_SN = SN_data["SN_observations"][ind]["sncosmo_model"].flux(current_date, IFS_args["waves"])
        else:
            f_lamb_SN = ones(len(IFS_args["waves"]), dtype=float64)*1e-22
        assert all(1 - isnan(f_lamb_SN))

        #args = cp.deepcopy(IFS_args)
        IFS_args["mdl"] = f_lamb_SN

        ETC_result = get_spec_with_err(redshift = 0, exp_time = row_to_add["exptime"], gal_flamb = SN_data["SN_observations"][ind]["gal_background"], bad_pixel_rate = survey_parameters["bad_pixel_rate"],
                                       show_plots = 0, **IFS_args)
        
        
        SN_data["SN_observations"][ind]["IFS_dates"].append(current_date)
        SN_data["SN_observations"][ind]["IFS_exptimes"].append(row_to_add["exptime"])

        noise = ETC_result["f_lamb_SN"]/ETC_result["spec_S/N"]
        assert all(1 - isnan(noise) - isinf(noise))
        SN_data["SN_observations"][ind]["IFS_fluxes"].append(
            ETC_result["f_lamb_SN"] + random.normal(size = len(ETC_result["f_lamb_SN"]))*noise)
        SN_data["SN_observations"][ind]["IFS_dfluxes"].append(noise)

    return SN_data

def run_observation_through_ground_ETC(SN_data, row_to_add, current_date, ground_depths):
    phases = (current_date - SN_data["SN_table"]["daymaxes"])/(1. + SN_data["SN_table"]["redshifts"])

    active_SN_mask = (phases >= -15.)*(phases <= 45.)
    
    inds = where(active_SN_mask)[0]
    
    for ind in inds:
        f_lamb_SN = SN_data["SN_observations"][ind]["sncosmo_model"].flux(current_date, ground_args["waves"])
        
        AB_mag = get_ground_AB_mag(f_lamb_SN, row_to_add["filt"])        
        
        AB_flux = 10.**(-0.4*(AB_mag - master_zp))
        AB_err = 0.2*10.**(-0.4*(ground_depths[row_to_add["filt"]] - master_zp))
        SNR_total = AB_flux/AB_err

    
        flux_with_noise = AB_flux + random.normal()*AB_flux/SNR_total
        SN_data["SN_observations"][ind]["dates"] = append(SN_data["SN_observations"][ind]["dates"], current_date)
        SN_data["SN_observations"][ind]["fluxes"] = append(SN_data["SN_observations"][ind]["fluxes"], flux_with_noise)
        SN_data["SN_observations"][ind]["dfluxes"] = append(SN_data["SN_observations"][ind]["dfluxes"], AB_flux/SNR_total)
        SN_data["SN_observations"][ind]["filts"] = append(SN_data["SN_observations"][ind]["filts"], row_to_add["filt"])
        SN_data["SN_observations"][ind]["ispar"] = append(SN_data["SN_observations"][ind]["ispar"], 0)

        
    return SN_data


def quantize_time(t):
    #return ceil(t/5.65001)*5.65 # The ...01 is for numerical accuracy.
    return ceil(array(t)/3.04001)*3.04 # The ...01 is for numerical accuracy.


def get_single_slew_or_filt_change(RAs, Decs, roll_angles, i, j, filt_inds, NN_filt_change, n_filters):
    filter_changes = filt_inds[j] - filt_inds[i]
    if filter_changes < 0:
        filter_changes += n_filters

    roll_angle_diff = abs(roll_angles[j] - roll_angles[i])*(180./pi)
    if roll_angle_diff > 180.:
        roll_angle_diff = 360. - roll_angle_diff

    return quantize_time(max(
        slew_fn(sqrt((RAs[i] - RAs[j])**2. + (Decs[i] - Decs[j])**2.)),
        slew_fn(roll_angle_diff),
        filter_changes*NN_filt_change
    ))


def get_slew_time(RAs, Decs, roll_angles, filt_inds = None, NN_filt_change = 120/7., show_solution = False, label_slews = True, label_filts = True,
                  square_degrees = 0, n_filters = 8, filt_names = None):
    if len(RAs) < 2:
        print("No slew")
        return 0, []
    elif len(RAs) == 2:
        print("Only one slew")
        dist = sqrt((RAs[1] - RAs[0])**2. + (Decs[1] - Decs[0])**2.)
        return quantize_time(slew_fn(dist)), [[RAs[0], Decs[0], RAs[1], Decs[1],
                                              dist,
                                              quantize_time(slew_fn(dist))]]

    if len(unique(roll_angles)) > 2:
        print("Large number of roll angles found!!! ", len(unique(roll_angles)))

    #show_solution = (len(RAs) > 20)
    if filt_inds == None:
        filt_inds = [0]*len(RAs)
    if filt_names == None:
        filt_names = ["None"]*len(RAs)

    valid_tmp_file = False
    attempts = 0

    while not valid_tmp_file:
        try:
            temp = open("tmp", 'w') #tempfile.NamedTemporaryFile()
            valid_tmp_file = True
        except:
            print("Couldn't make temp file!")
            attempts += 1
            assert attempts < 100
            time.sleep(10)

    print("temp", temp)

    temp.write("TYPE: TSP\n")
    temp.write("DIMENSION: %i\n" % (len(RAs) + 1))
    temp.write("EDGE_WEIGHT_TYPE: EXPLICIT\n")
    temp.write("EDGE_WEIGHT_FORMAT: UPPER_ROW\n")
    temp.write("EDGE_WEIGHT_SECTION\n")
    
    for i in range(len(RAs)):
        for j in range(i+1,len(RAs)):
            temp.write("%0.f " % (get_single_slew_or_filt_change(RAs, Decs, roll_angles, i, j, filt_inds, NN_filt_change, n_filters)*100.))
        temp.write('0 \n')
    temp.write("EOF\n")

    temp.flush()

    cmd = "/Users/" + whoami + "/Dropbox/Shared/concorde/TSP/concorde -x tmp "# + temp.name # -x deletes tmp files
    print(cmd)
    print(subprocess.getoutput(cmd))
    
    
    f = open(temp.name.split("/")[-1] + ".sol", 'r')
    lines = f.read().split(None)
    f.close()
    subprocess.getoutput("rm -f " + temp.name.split("/")[-1] + ".sol")
    temp.close()
    

    assert int(lines[0]) == len(RAs) + 1, str(lines)
    del lines[0]

    sol = [int(item) for item in lines]
    ind = sol.index(len(sol) - 1)
    sol = sol[ind+1:] + sol[:ind]

    labeled = []
    total_time = 0.


    if show_solution:
        if label_filts:
            for i in range(len(RAs)):
                plt.plot(RAs[i], Decs[i], '.',
                         color = {"R087": (1, 0.5, 1), "Z087": 'm', "Y106": 'b', "J129": 'c', "H158": 'g', "F184": 'r', "K213": 'k', "W149": 'gray', "W146": 'gray', "None": 'gray', "Pointings": 'gray'}[filt_names[i]],
                         label = filt_names[i]*(labeled.count(filt_names[i]) == 0)*label_filts, zorder = 3)
                labeled.append(filt_names[i])
        else:
            plt.plot(RAs, Decs, '.', color = 'b')


        circle1=plt.Circle((0,0),sqrt(square_degrees/pi), fill = False, label = 'Field Edge')
        fig = plt.gcf()
        fig.gca().add_artist(circle1)

        patches = []
        for i in range(len(RAs)):
            xs, ys = field_outline(RAs_to_add[i], Decs_to_add[i], WFI_orient = 0)
            
            """
            xs = concatenate((xs, xs[:,:1]), axis = 1)
            ys = concatenate((ys, ys[:,:1]), axis = 1)
            
            for j in range(len(xs)):
                plt.plot(xs[j], ys[j], linewidth = 0.25, color = [(0.75, 0.75, 0.75), (0.25, 0.25, 0.25), (0.5,0.5,0.5)][i%3])
            """
            for j in range(len(xs)):
                poly = Polygon(list(zip(xs[j],ys[j])),facecolor=[(1, 0.8, 0.8), (0.8, 1, 0.8), (0.8, 0.8, 1.)][i%3],edgecolor='none')
                plt.gca().add_patch(poly)

    print("slew ", end=' ') 


    ordered_sol = []

    slew_times_found = []
    for i in range(len(RAs) - 1):
        slew_time = get_single_slew_or_filt_change(RAs, Decs, roll_angles, sol[i], sol[i+1], filt_inds, NN_filt_change, n_filters)
        slew_times_found.append(slew_time)
    slew_times_found = sort(list(set(slew_times_found)))

    color_dict = {}

    if label_slews and show_solution:
        for i in range(len(slew_times_found)):
            if len(slew_times_found) <= 4:
                color_dict["%.2f" % slew_times_found[i]] = ['b', 'g', 'orange', 'r'][i]
            else:
                color_dict["%.2f" % slew_times_found[i]] = ['m', 'b', 'c', 'g', 'orange', 'r', 'k'][i]

    for i in range(len(RAs) - 1):
        slew_time = get_single_slew_or_filt_change(RAs, Decs, roll_angles, sol[i], sol[i+1], filt_inds, NN_filt_change, n_filters)
        ordered_sol.append([RAs[sol[i]], Decs[sol[i]], RAs[sol[i+1]], Decs[sol[i+1]],
                            sqrt((RAs[sol[i]] - RAs[sol[i+1]])**2. + (Decs[sol[i]] - Decs[sol[i+1]])**2.),
                            slew_time])

        print(slew_time, end=' ') 
        total_time += slew_time
        
        if show_solution:
            txt_slew = "%.2f" % slew_time

            plt.plot([RAs[sol[i]], RAs[sol[i+1]]],
                     [Decs[sol[i]], Decs[sol[i+1]]],
                     color = color_dict.get(txt_slew, 'k'), label = (txt_slew + "s")*(labeled.count(txt_slew) == 0)*label_slews)
            labeled.append(txt_slew)
            #plt.text(mean([RAs[sol[i]], RAs[sol[i+1]]]), mean([Decs[sol[i]], Decs[sol[i+1]]]), "    %.2f" % slew_time, size = 6)


    if show_solution:
        plt.legend(loc = 'best', fontsize = 10)
        plt.axes().set_aspect('equal', 'datalim')
        plt.axes().set_aspect('equal', 'datalim')
        plt.axis('off')
        
        #plt.title("t=%.3f" % total_time)
        plt.savefig("slew_solution.eps", bbox_inches = 'tight')
        plt.close()
        
    print("done ", total_time)
    
    return total_time, ordered_sol


"""

slew_fn = file_to_fn(wfirst_data_path + "/pixel-level/input/slewtable_I32000_t0.40_vmax0.12_Tr2.0.txt", col = 2)
rs = sqrt(random.random(size = 20))*(sqrt(5./pi) - sqrt(0.28)/2)
ths = random.random(size = len(rs))*2*pi
xs = rs*cos(ths)
ys = rs*sin(ths)
roll_angles = random.randint(2, size = len(rs))*5.*(pi/180.)

filt_inds = [0]*len(xs)
get_slew_time(xs, ys, roll_angles, filt_inds = filt_inds, show_solution = True, label_filts = True, label_slews = False, square_degrees = 5)


if 0:
    filt_inds = [1]*17 + [2]*17 + [3]*17
    filt_names = ["Z087"]*17 + ["H158"]*17 + ["F184"]*17
else:
    filt_inds = [1]*10 + [3]*10
    filt_names = ["H158"]*10 + ["F184"]*10
get_slew_time(xs, ys, roll_angles, filt_inds = filt_inds, show_solution = True, label_filts = True, label_slews = False, filt_names = filt_names, square_degrees = 5)
stophere
"""


def add_observations_to_sndata(SN_data, rows_to_add, current_date, ground_depths, ordered_sol):
    rows_added = Table(names = ["date", "filt", "exptime", "RA", "dec", "orient", "instr", "SNind"],
                       dtype= ("f8", "S10", "f8", "f8", "f8", "f8", "S10", "i4"))

    filt_set = list(set(list(rows_to_add["filt"])))
    time_used = 0.
    slew_time = 0.

    for filt in filt_set:
        inds = where(
            (rows_to_add["filt"] == filt)*(abs(rows_to_add["date"] - current_date) < 1.)
                     )[0]
        print(filt, inds)
        inds = sort(inds)[::-1]
        print(inds)
        if len(inds) > 0 and rows_to_add["instr"][inds[0]] != "ground":
            tmp_slew_time, tmp_ordered_sol = get_slew_time(rows_to_add["RA"][inds], rows_to_add["dec"][inds], rows_to_add["orient"][inds])
            tmp_slew_time += quantize_time(3. * 120./7.) # 51 seconds for changing filters (three changes)

            for item in tmp_ordered_sol:
                ordered_sol.add_row([current_date, filt] + item)

            time_used += tmp_slew_time
            slew_time += tmp_slew_time
        
        for ind in inds:
            # For this filter, for current_date, for row "ind"
            if rows_to_add["instr"][ind] != "ground":
                time_used += quantize_time(rows_to_add["exptime"][ind])
                SN_data = run_observation_through_ETC(SN_data, rows_to_add[ind], current_date)
            else:
                SN_data = run_observation_through_ground_ETC(SN_data, rows_to_add[ind], current_date, ground_depths)

            rows_added.add_row(rows_to_add[ind])
            rows_to_add.remove_row(ind) # Done with this observation!

        # Finished adding all data for this filter, now, take variance-weighted mean of all repeat observations for a given SN (either because of dithers, or parallels)
        SN_data = compress_lc_data(SN_data, current_date, filt)

    return SN_data, rows_to_add, rows_added, time_used, slew_time, ordered_sol


def find_SNe(SN_data, current_date):
    for i in range(SN_data["nsne"]):
        if SN_data["SN_observations"][i]["found_date"] == None:
            if len(SN_data["SN_observations"][i]["fluxes"]) > 0:
                inds = where(abs(SN_data["SN_observations"][i]["dates"] - current_date) < 1)
                SNRs = SN_data["SN_observations"][i]["fluxes"][inds]/SN_data["SN_observations"][i]["dfluxes"][inds]
                SNRs = sort(SNRs)
                if len(SNRs) > 1 and SNRs[-2] >= 4:
                    SN_data["SN_observations"][i]["found_date"] = current_date
    return SN_data


def get_IFS_exptimes(redshift, sncosmo_model, daymax, IFS_trigger_params, gal_flamb):
    IFS_exptimes = {}

    for SNR in IFS_trigger_params["SNR_set"]:
        IFS_args["restframe_bins"] = [IFS_trigger_params["normalization_wavelength_range"][0]*(1. + redshift), IFS_trigger_params["normalization_wavelength_range"][1]*(1 + redshift)]

        flamb_model = sncosmo_model.flux(daymax, IFS_args["waves"])
        IFS_args["mdl"] = flamb_model

        IFS_exptimes[SNR] = quantize_time(solve_for_exptime(S_to_N = SNR, redshift = 0, gal_flamb = gal_flamb, bad_pixel_rate = 0,
                                                            key1 = "rest_frame_band_S/N",
                                                            key2 = (IFS_trigger_params["normalization_wavelength_range"][0]*(1. + redshift),
                                                                    IFS_trigger_params["normalization_wavelength_range"][1]*(1 + redshift)), **IFS_args)
                                          )

    print("Found ", redshift, IFS_exptimes)
    return IFS_exptimes










def plan_and_add_triggers(SN_data, rows_to_add, current_date, cadence, IFS_trigger_params, parallel_filters, square_degrees):
    """This function spreads parallel observations over the sky."""


    trigger_prob_fn = interp1d(IFS_trigger_params["trigger_redshifts"],
                               clip([item*IFS_trigger_params["trigger_scaling"] for item in IFS_trigger_params["trigger_fraction"]], 0, 1),
                               kind = 'linear')


    possible_inds = array([], dtype=int32)
    possible_metrics = array([], dtype=int32) # In this case, linear degrees from other visits. Greedy algorithm, but better than nothing.
    possible_redshifts = array([], dtype=float64)

    for i in range(SN_data["nsne"]):
        if SN_data["SN_observations"][i]["found_date"] == current_date: # Newly found SN, candidate for trigger

            possible_inds = append(possible_inds, i)
            possible_redshifts = append(possible_redshifts, SN_data["SN_table"]["redshifts"][i])


            WFI_RA, WFI_dec = IFS_to_WFI(SN_data["SN_table"]["RAs"][i], SN_data["SN_table"]["Decs"][i], date_to_roll_angle(SN_data["SN_table"]["daymaxes"][i]))

            if WFI_RA**2. + WFI_dec**2. > square_degrees/pi:
                possible_metrics = append(possible_metrics, random.random() - 2)
            else:
                possible_metrics = append(possible_metrics, random.random())
    
    
    # Now, randomly choose SNe to trigger
    print("current_date, possible_metrics ", current_date, possible_metrics)


    for i in range(len(possible_metrics)):
        quantile = percentileofscore(possible_metrics, possible_metrics[i])/100.
        select_this_SN = 1 - quantile < trigger_prob_fn(possible_redshifts[i]) # Has to be less than, or else 0 is always triggered on!
        print(possible_redshifts[i], possible_metrics[i], "selected ", select_this_SN)

        if select_this_SN:
            if IFS_trigger_params["adjust_each_SN_exp_time"]:
                IFS_exptimes = get_IFS_exptimes(redshift = possible_redshifts[i],
                                                sncosmo_model = SN_data["SN_observations"][possible_inds[i]]["sncosmo_model"], gal_flamb = SN_data["SN_observations"][possible_inds[i]]["gal_background"],
                                                daymax = SN_data["SN_table"]["daymaxes"][possible_inds[i]], IFS_trigger_params = IFS_trigger_params)
                
            else:
                IFS_exptimes = IFS_trigger_params["uniform_IFS_exptimes"][possible_redshifts[i]]

            
            print("IFS_exptimes ", IFS_exptimes)
                
            trigger_dates = []

            for phase in IFS_trigger_params["phases"]:
                if trigger_dates == []:
                    next_date = max(current_date + cadence, SN_data["SN_table"]["daymaxes"][possible_inds[i]] + phase*(1. + possible_redshifts[i]))
                else:
                    next_date = max(trigger_dates[-1] + cadence, SN_data["SN_table"]["daymaxes"][possible_inds[i]] + phase*(1. + possible_redshifts[i]))
                trigger_dates.append(next_date)

            # Add references
            trigger_SNRs = IFS_trigger_params["SNRs"] + [IFS_trigger_params["SNR_set"][-1]]*IFS_trigger_params["number_of_reference_dithers"]

            dates_wrt_one_year = linspace(0., 365./2., IFS_trigger_params["number_of_reference_dithers"] + 1.)
            dates_wrt_one_year = dates_wrt_one_year[:-1]
            dates_wrt_one_year -= median(dates_wrt_one_year)
            print("dates_wrt_one_year", dates_wrt_one_year)

            for dwrtoy in dates_wrt_one_year:
                trigger_dates.append(SN_data["SN_table"]["daymaxes"][possible_inds[i]] + 365.*(1. + possible_redshifts[i]) + dwrtoy)
                
            trigger_dates = round_to_cadence(array(trigger_dates), cadence)
            print("trigger_dates, trigger_SNRs", trigger_dates, trigger_SNRs)
            

            for date, SNR in zip(trigger_dates, trigger_SNRs):
                WFI_RA, WFI_dec = IFS_to_WFI(SN_data["SN_table"]["RAs"][possible_inds[i]], SN_data["SN_table"]["Decs"][possible_inds[i]], date_to_roll_angle(date))
                rows_to_add.add_row((date, "W149", IFS_exptimes[SNR], WFI_RA, WFI_dec, date_to_roll_angle(date), "WFI", possible_inds[i])) # If you see W149 in a LC, something is wrong!



            r = sqrt(random.random())*sqrt(square_degrees/pi)
            theta = random.random()*2*pi

            CC_RA = r*cos(theta)
            CC_Dec = r*sin(theta)

            for date, phase, SNR in zip(trigger_dates[:2], IFS_trigger_params["phases"], trigger_SNRs[:2]):
                # For every SN Ia triggered on, assume there is one CC SN that is falsely triggered on for two epochs. The assumption of a random place is conservative; we're not doing anything with the parallels.
                CC_WFI_RA, CC_WFI_Dec = IFS_to_WFI(CC_RA, CC_Dec, date_to_roll_angle(date))

                if phase < 3:
                    if random.random() < 0.2:
                        rows_to_add.add_row((date, "W149", IFS_exptimes[SNR], CC_WFI_RA, CC_WFI_Dec, date_to_roll_angle(date), "WFI", -2))
                        print("Triggering on CC SN, phase=", phase)
            


    return rows_to_add


def optimize_triggers(rows_to_add, current_date, cadence, parallel_filters):

    #     rows_to_add = Table(names = ["date", "filt", "exptime", "RA", "dec", "orient", "instr", "SNind"],
    inds = where((rows_to_add["date"] == current_date + cadence)*(rows_to_add["filt"] == "W149")) # Plan the next visit
    if len(inds[0]) > 0:
        print(inds)
        print(rows_to_add[inds])
        RA_need_filt = rows_to_add[inds]["RA"]
        Dec_need_filt = rows_to_add[inds]["dec"]
        exp_need_filt = rows_to_add[inds]["exptime"]

        print(RA_need_filt)
        print(Dec_need_filt)

        best_metric = -1
        best_filts = None
        for i in range(400):
            possible_choice = random.choice(parallel_filters, size=len(RA_need_filt), replace=True)
            #print possible_choice
            this_metric = 1e10
            
            for filt in unique(possible_choice):
                filt_inds = where(possible_choice == filt)

                RA_this_filt = RA_need_filt[filt_inds]
                Dec_this_filt = Dec_need_filt[filt_inds]
                
                for j in range(len(RA_this_filt)):
                    for k in range(j+1, len(RA_this_filt)):
                        this_metric = min(this_metric, (RA_this_filt[j] - RA_this_filt[k])**2. + (Dec_this_filt[j] - Dec_this_filt[k])**2.
                                          )
                        #print this_metric

            for filt in parallel_filters:
                actual_frac = sum(possible_choice == filt)/float(len(possible_choice))
                actual_frac = clip(actual_frac, 0.001, 1.) # Prevent zeros

                target_frac = 1./len(parallel_filters)

                this_metric -= abs(log(target_frac/actual_frac))

            if parallel_filters.count("F184") and len(parallel_filters) > 1:
                target_for_F184 = 3./(len(parallel_filters) + 2.)
                actual_frac = sum(exp_need_filt[where(possible_choice == "F184")])/sum(exp_need_filt)

                actual_frac = clip(actual_frac, 0.001, 1.) # Prevent zeros

                this_metric -= 2.*abs(log(target_for_F184/actual_frac))

            if this_metric > best_metric:
                best_filts = possible_choice
                best_metric = this_metric
                print("New best", best_filts, best_metric)

        for i, ind in enumerate(inds[0]):
            rows_to_add[ind]["filt"] = best_filts[i]
                


    return rows_to_add



def plan_and_add_triggers_old(SN_data, rows_to_add, current_date, cadence, IFS_trigger_params, parallel_filters, square_degrees):
    """This function targets SNe in the parallels, but rolling quickly doesn't work. So the assumption of a constant roll angle is off. I'm replacing this with a function that just spreads observations spatially."""


    trigger_roll_angle = date_to_roll_angle(current_date + cadence + 15.) # can hold a roll angle for 30 days

    active_SNe = array([SN_data["SN_observations"][i]["found_date"] != None for i in range(SN_data["nsne"])]) # SNe that have been discovered, and are no later than 10 rest-frame days after maximum.
    active_SNe *= (current_date - SN_data["SN_table"]["daymaxes"])/(1. + SN_data["SN_table"]["redshifts"]) <= 10.
    active_inds = where(active_SNe)


    possible_inds = array([], dtype=int32)
    possible_NSNe = array([], dtype=int32)
    possible_redshifts = array([], dtype=float64)

    for i in range(SN_data["nsne"]):
        if SN_data["SN_observations"][i]["found_date"] == current_date: # Newly found SN, candidate for trigger

            possible_inds = append(possible_inds, i)
            possible_redshifts = append(possible_redshifts, SN_data["SN_table"]["redshifts"][i])


            WFI_RA, WFI_dec = IFS_to_WFI(SN_data["SN_table"]["RAs"][i], SN_data["SN_table"]["Decs"][i], trigger_roll_angle)

            if sum(active_SNe) > 0:
                inside_field_mask = inside_field(SN_data["SN_table"]["RAs"][active_inds], SN_data["SN_table"]["Decs"][active_inds], WFI_RA, WFI_dec, trigger_roll_angle)
                possible_NSNe = append(possible_NSNe, sum(inside_field_mask))
            else:
                possible_NSNe = append(possible_NSNe, 0)

            """
            if possible_NSNe[-1] > 4:
                inside_field_inds = where(inside_field_mask)
                plt.plot(SN_data["SN_table"]["RAs"][i], SN_data["SN_table"]["Decs"][i], 'o', color = 'b')
                plt.plot(SN_data["SN_table"]["RAs"][active_inds][inside_field_inds], SN_data["SN_table"]["Decs"][active_inds][inside_field_inds], '.', color = 'g')
                
                outside_field_inds = where(1 - inside_field_mask)
                plt.plot(SN_data["SN_table"]["RAs"][active_inds][outside_field_inds], SN_data["SN_table"]["Decs"][active_inds][outside_field_inds], '.', color = 'r')
                plt.title("Date %.1f IFS RA %.3f Dec %.3f WFI RA %.3f Dec %.3f" % (current_date, SN_data["SN_table"]["RAs"][i], SN_data["SN_table"]["Decs"][i], WFI_RA, WFI_dec))
                plt.axes().set_aspect('equal', 'datalim')

                plt.savefig("parallel_test.pdf")
                plt.close()
                stophere"""
    
    
    # Now, randomly choose SNe to trigger
    print("current_date, possible_NSNe ", current_date, possible_NSNe)


    for i in range(len(possible_NSNe)):
        quantile = percentileofscore(possible_NSNe, possible_NSNe[i])/100.
        #select_this_SN = 1 - quantile < trigger_prob_fn(possible_redshifts[i]) # Has to be less than, or else 0 is always triggered on!
        print(possible_redshifts[i], possible_NSNe[i], "selected ", select_this_SN)

        if select_this_SN:
            if IFS_trigger_params["adjust_each_SN_exp_time"]:
                IFS_exptimes = get_IFS_exptimes(redshift = possible_redshifts[i],
                                                sncosmo_model = SN_data["SN_observations"][possible_inds[i]]["sncosmo_model"], gal_flamb = SN_data["SN_observations"][possible_inds[i]]["gal_background"],
                                                daymax = SN_data["SN_table"]["daymaxes"][possible_inds[i]], IFS_trigger_params = IFS_trigger_params)
                
            else:
                IFS_exptimes = IFS_trigger_params["uniform_IFS_exptimes"][possible_redshifts[i]]

            
            print("IFS_exptimes ", IFS_exptimes)
                
            trigger_dates = []

            for phase in IFS_trigger_params["phases"]:
                if trigger_dates == []:
                    next_date = max(current_date + cadence, SN_data["SN_table"]["daymaxes"][possible_inds[i]] + phase*(1. + possible_redshifts[i]))
                else:
                    next_date = max(trigger_dates[-1] + cadence, SN_data["SN_table"]["daymaxes"][possible_inds[i]] + phase*(1. + possible_redshifts[i]))
                trigger_dates.append(next_date)

            trigger_SNRs = IFS_trigger_params["SNRs"] + [IFS_trigger_params["SNR_set"][-1]]
            trigger_dates.append(SN_data["SN_table"]["daymaxes"][possible_inds[i]] + 365.)
            
            trigger_dates = round_to_cadence(array(trigger_dates), cadence)
            print("trigger_dates, trigger_SNRs", trigger_dates, trigger_SNRs)

            WFI_RA, WFI_dec = IFS_to_WFI(SN_data["SN_table"]["RAs"][possible_inds[i]], SN_data["SN_table"]["Decs"][possible_inds[i]], trigger_roll_angle)

            for date, SNR in zip(trigger_dates, trigger_SNRs):
                rows_to_add.add_row((date, random.choice(parallel_filters), IFS_exptimes[SNR], WFI_RA, WFI_dec, trigger_roll_angle, "WFI", possible_inds[i]))



            r = sqrt(random.random())*sqrt(square_degrees/pi)
            theta = random.random()*2*pi

            CC_RA = r*cos(theta)
            CC_Dec = r*sin(theta)

            for date, phase, SNR in zip(trigger_dates, IFS_trigger_params["phases"], trigger_SNRs):
                # For every SN Ia triggered on, assume there is one CC SN that is falsely triggered on for two epochs. The assumption of a random place is conservative; we're not doing anything with the parallels.
                if phase < 0:
                    rows_to_add.add_row((date, random.choice(parallel_filters), IFS_exptimes[SNR], CC_RA, CC_Dec, trigger_roll_angle, "WFI", -2))
                    print("Triggering on CC SN, phase=", phase)
            
    
    return rows_to_add
    

def get_RA_dec_for_cadence(square_degrees, next_angle):
    #small_step = sqrt(0.28*0.5)
    #large_step = sqrt(0.28*2)
    small_step = 0.37628
    large_step = 0.82034
    # Fill factor is about 90%

    field_radius = sqrt(square_degrees/pi)

    x_vals = arange(around(2*field_radius/small_step))*small_step
    x_vals -= median(x_vals)
    #y_vals = arange(-around(field_radius*2.5/large_step), around(field_radius*2.5/large_step) + 0.001)*large_step + small_step

    


    RAs_to_add = []
    Decs_to_add = []

    for xv in x_vals:
        if abs(xv) <= field_radius:
            ylength = sqrt(field_radius**2. - xv**2.)*2.
            rounded_size = around(ylength/large_step + 0.1) # 0.1 to give the benefit of the doubt and minimize edge effects.
            
            """
            # I thought this would be a good idea, but we lose a lot at the edges of the field
            if around(2*field_radius/large_step) % 2 == 0:
                rounded_size = 2.*around(rounded_size/2.)
            else:
                rounded_size = 2.*around((rounded_size - 1.)/2.) + 1
            """

            y_vals = arange(rounded_size)*large_step
            y_vals -= median(y_vals)

            for yv in y_vals:
                x0 = xv
                y0 = yv
                
                RAs_to_add.append(cos(next_angle)*x0 - y0*sin(next_angle))
                Decs_to_add.append(sin(next_angle)*x0 + y0*cos(next_angle))
    return RAs_to_add, Decs_to_add

"""
square_degrees = 9.6

RAs_to_add, Decs_to_add = get_RA_dec_for_cadence(square_degrees = square_degrees, next_angle = 0)

SN_rs = sqrt(random.random(size = 350))*sqrt(square_degrees/pi)
SN_ths = random.random(size = 350)*2*pi

SN_RAs = SN_rs*cos(SN_ths)
SN_decs = SN_rs*sin(SN_ths)


dRAs = arange(-0.2, 0.2001, 0.01)
ddecs = arange(-0.2, 0.2001, 0.01)

number_im = zeros([len(dRAs), len(ddecs)])

for j in range(len(dRAs)):
    print j, len(dRAs)

    for k in range(len(ddecs)):
        mask = 0
        for i in range(len(RAs_to_add)):
            mask += inside_field(SN_RAs, SN_decs, RAs_to_add[i] + dRAs[j], Decs_to_add[i] + ddecs[k], WFI_orient=0)
    
        mask = mask > 0
        print sum(mask)
        
        number_im[j,k] = sum(mask)


plt.imshow(number_im, interpolation = 'nearest')
plt.colorbar()
plt.savefig("tmp.pdf")
plt.close()

        
stop_here
"""

"""
for square_degrees in arange(1., 20.1, 1.):
    print square_degrees
    RAs_to_add, Decs_to_add = get_RA_dec_for_cadence(square_degrees = square_degrees, next_angle = 0)
    
    circle1=plt.Circle((0,0),sqrt(square_degrees/pi), fill = False)
    fig = plt.gcf()
    fig.gca().add_artist(circle1)
    for i in range(len(RAs_to_add)):
        xs, ys = field_outline(RAs_to_add[i], Decs_to_add[i], WFI_orient = 0)

        xs = concatenate((xs, xs[:,:1]), axis = 1)
        ys = concatenate((ys, ys[:,:1]), axis = 1)

        for j in range(len(xs)):
            plt.plot(xs[j], ys[j], linewidth = 0.25, color = 'rgbkcm'[i%6])


    plt.xlim(-5, 5)
    plt.ylim(-5, 5)

    #plt.plot(RAs_to_add, Decs_to_add, '.')
    plt.axes().set_aspect('equal', 'datalim')

    plt.savefig("square_degrees=%04i.pdf" % (square_degrees*100))
    plt.close()

    slew_fn = file_to_fn(wfirst_data_path + "/pixel-level/input/slew_settle_170419.txt", col = 2)
    get_slew_time(RAs_to_add, Decs_to_add, roll_angles = [0]*len(RAs_to_add), show_solution = 1, square_degrees = square_degrees, filt_names = ["Pointings"]*len(RAs_to_add))
    commands.getoutput("mv slew_solution.eps slew_solution_%04i.eps" % (square_degrees*100))


stop_here
"""


def plan_and_add_cadence(SN_data, wfi_time_left, current_date, rows_to_add,
                         square_degrees, cadence,
                         tier_filters, tier_exptimes, tier_dithers):

    next_angle = date_to_roll_angle(current_date + cadence)
    wide_gap_degrees = 8.546/40.88 * (0.11*4096/3600.) # 8.546 mm

    dither_positions = []
    for i in arange(-3., 4):
        for j in arange(-3., 4):
            dither_positions.append([i**2. + j**2., i*wide_gap_degrees, j*wide_gap_degrees])

    dither_positions.sort()


    current_count = 0

    for filt, expt, dith in zip(tier_filters, tier_exptimes, tier_dithers):
        if len(filt) > 1:
            # If WFI, not ground
            if wfi_time_left:
                RAs_to_add, Decs_to_add = get_RA_dec_for_cadence(square_degrees, next_angle)
                
                for i in range(len(RAs_to_add)):
                    for j in range(dith):
                        rows_to_add.add_row((current_date + cadence, filt, expt,
                                             RAs_to_add[i] + dither_positions[current_count][1],
                                             Decs_to_add[i] + dither_positions[current_count][2], next_angle, "WFI", -1))

                current_count += 1

            else:
                pass # Out of time, nothing to add
        else:
            rows_to_add.add_row((current_date + cadence, filt, 0, 0, 0, 0, "ground", -1))

    return rows_to_add



def plan_and_add_WFI_cadence_rectangle(SN_data, current_date, rows_to_add,
                                       RA_linear_degrees, Dec_linear_degrees,
                                       tier_filters, tier_exptimes, dithers):

    small_step = sqrt(0.28*0.5)
    large_step = sqrt(0.28*2)


    orient = (around(sun_angle/(pi/2.)) % 4)*(pi/2.)

    if abs(orient) < 0.01 or abs(orient - pi) < 0.01:
        x_vals = arange(around(RA_linear_degrees/small_step))*small_step
        y_vals = arange(around(Dec_linear_degrees/large_step))*large_step
    elif abs(orient - pi/2.) < 0.01 or abs(orient - 1.5*pi) < 0.01:
        x_vals = arange(around(RA_linear_degrees/large_step))*large_step
        y_vals = arange(around(Dec_linear_degrees/small_step))*small_step
    else:
        raise Exception("Weird orient!", orient)


    wide_gap_degrees = 8.546/40.88 * (0.11*4096/3600.) # 8.546 mm

    for xv in x_vals:
        for yv in y_vals:
            for filt, expt in zip(tier_filters, tier_exptimes):
                for i in range(dithers):
                    rows_to_add.add_row((current_date + cadence, filt, expt, xv + i*wide_gap_degrees, yv + i*wide_gap_degrees, orient, "WFI", -1))

    return rows_to_add


def run_survey(SN_data, square_degrees, tier_filters, tier_exptimes, ground_depths, dithers, cadence, fraction_time, total_survey_time, hours_per_visit, IFS_trigger_params, parallel_filters, survey_duration):

    starting_time = 31557000*total_survey_time*fraction_time  # Seconds in total_survey_time years
    total_time_left = starting_time
    total_slew_time = 0.0

    # Start with empty table
    # instr is WFI or ground
    # SNind is index of SN (if any) in IFS
    # RA, dec are WFI, not IFS

    observation_table = Table(names = ["date", "filt", "exptime", "RA", "dec", "orient", "instr", "SNind"],
                              dtype= ("f8", "S10", "f8", "f8", "f8", "f8", "S10", "i4"))
    ordered_sol = Table(names = ["date", "filt", "RA1", "dec1", "RA2", "dec2", "deg", "time"],
                        dtype= ("f8", "S10", "f8", "f8", "f8", "f8", "f8", "f8"))

    # RA is x, Dec is y

    print(observation_table)
    current_date = 0.
    current_trigger_scaling = 1.

    # Start with empty table
    rows_to_add = Table(names = ["date", "filt", "exptime", "RA", "dec", "orient", "instr", "SNind"],
                        dtype= ("f8", "S10", "f8", "f8", "f8", "f8", "S10", "i4"))

    SN_data["time_remaining_values"] = [[(current_date, total_time_left, total_time_left, current_trigger_scaling)]]

    while (len(rows_to_add["filt"]) > 0 or current_date == 0) and current_date < 365*(3 + survey_duration):
        # The extra 1.5 on the 2 years is to get references + have some margin
        IFS_trigger_params["trigger_scaling"] = current_trigger_scaling


        # Step 1: Add previously planned observations
        SN_data, rows_to_add, rows_added, time_used, slew_time, ordered_sol = add_observations_to_sndata(SN_data, rows_to_add, current_date,
                                                                                                         ground_depths, ordered_sol = ordered_sol)
        for i in range(len(rows_added)):
            observation_table.add_row(rows_added[i])
        total_time_left -= time_used
        total_slew_time += slew_time
        
        # Step 2: Find SNe in observations
        SN_data = find_SNe(SN_data, current_date)
        

        # Step 3: Plan triggers
        if total_time_left > sum(rows_to_add["exptime"] + 50.):
            rows_to_add = plan_and_add_triggers(SN_data, rows_to_add, current_date, cadence, IFS_trigger_params = IFS_trigger_params,
                                                parallel_filters = parallel_filters, square_degrees = square_degrees)
        rows_to_add = optimize_triggers(rows_to_add, current_date, cadence, parallel_filters)

        # Step 4: Plan cadenced observations for next time

        estimated_time_left = total_time_left - sum(rows_to_add["exptime"] + 50.)
        wfi_time_left = estimated_time_left > 0
        # This is an approximate ending condition; we'll see how close total_time_left comes to zero.
            
        rows_to_add = plan_and_add_cadence(SN_data = SN_data, wfi_time_left = wfi_time_left, current_date = current_date,
                                           rows_to_add = rows_to_add,
                                           square_degrees = square_degrees, cadence = cadence,
                                           tier_filters = tier_filters, tier_exptimes = tier_exptimes, tier_dithers = dithers)
            
        
        
        print("\nEnd of day", current_date, "time left", total_time_left, time.asctime(), '\n')

        current_date += cadence
        SN_data["time_remaining_values"][0].append((current_date, total_time_left, estimated_time_left, current_trigger_scaling))
    
        """
        fraction_of_time_used = (starting_time - estimated_time_left)/starting_time
        fraction_of_dates_used = current_date/(survey_duration*365.)


        raw_rate = fraction_of_dates_used/(fraction_of_time_used + 0.00001) # Just to avoid divide by zero
        # Now, it's time to smooth a bit
        rate_scaling = raw_rate**1.0 # Actually, no smoothing works well
        rate_scaling = min(rate_scaling, 1./max(IFS_trigger_params["trigger_fraction"]))
        IFS_trigger_params["trigger_fraction"] = [item*rate_scaling for item in IFS_trigger_params["trigger_fraction"]]
        print "heading into, rate_scaling ", current_date, rate_scaling, raw_rate, IFS_trigger_params["trigger_fraction"], fraction_of_time_used, fraction_of_dates_used
        """

        if current_date > 10*cadence:
            # Let's take a look at the last five steps
            time_in_last_five = SN_data["time_remaining_values"][0][-6][1] - SN_data["time_remaining_values"][0][-1][1]
            relative_time_scaling = hours_per_visit*fraction_time*3600*5./(time_in_last_five + 1e-6)

            current_trigger_scaling = relative_time_scaling*mean([SN_data["time_remaining_values"][0][i][3] for i in range(-6, 0)])

            current_trigger_scaling = min(current_trigger_scaling, 1./min([item for item in IFS_trigger_params["trigger_fraction"] if item > 0]))

        print("date ", current_date, "current_trigger_scaling ", current_trigger_scaling, "total_time_left ", total_time_left)


    print("Done with tier. Slew time = ", total_slew_time)
    SN_data["total_time_used"] = starting_time - total_time_left
    SN_data["observation_table"] = observation_table
    SN_data["ordered_sol"] = ordered_sol
    return SN_data




def make_SNe(square_degrees, cadence, survey_duration, hours_per_visit, rates_fn, redshift_set, IFS_trigger_params, tier_filters, tier_exptimes, ground_depths, dithers, parallel_filters,
             fraction_time, total_survey_time, redshift_step = 0.05, salt2_model = True, verbose = False, phase_buffer = 20, survey_fields = "None"):
    assert square_degrees < 500, "Should use more accurate formula for large surveys!"


    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    frac_of_sky = (square_degrees)/(4.*pi*(180./pi)**2.)
    
    if verbose:
        print("frac_of_sky ", frac_of_sky)
        

    volume_in_shells = cosmo.comoving_volume(redshift_step/2. + redshift_set).value - cosmo.comoving_volume(redshift_set - redshift_step/2.).value
    SNe_in_survey_field = volume_in_shells*frac_of_sky*rates_fn(redshift_set) * survey_duration/(1. + redshift_set) * 1e-4

    if verbose:
        print(SNe_in_survey_field, sum(SNe_in_survey_field))

    SNe_actual = [random.poisson(item) for item in SNe_in_survey_field]
    if verbose:
        print(SNe_actual, sum(SNe_actual))
    redshifts = []
    for i in range(len(redshift_set)):
        redshifts += [redshift_set[i]]*SNe_actual[i]
    redshifts = array(redshifts)

    assert len(redshifts) == sum(SNe_actual)

    daymaxes = random.random(size = sum(SNe_actual))*survey_duration*365
    not_at_edge = where((daymaxes/(1. + redshifts) > phase_buffer)*((survey_duration*365 - daymaxes)/(1. + redshifts) > phase_buffer))

    redshifts = redshifts[not_at_edge]
    daymaxes = daymaxes[not_at_edge]
    nsne = len(redshifts)

    if verbose:
        print("SNe not near beginning/end: ", nsne)

    rs = sqrt(random.random(size = nsne))*sqrt(square_degrees/pi)
    thetas = random.random(size = nsne)*2*pi
        
    RAs = rs*cos(thetas)
    Decs = rs*sin(thetas)


    rs = sqrt(random.random(size = 10))*sqrt(square_degrees/pi)
    thetas = random.random(size = 10)*2*pi
        
    test_RAs = rs*cos(thetas)
    test_Decs = rs*sin(thetas)
    test_Mags = arange(24., 31.1, 0.4)

    test_points = {"RAs": [], "Decs": [], "Mags": [], "Observations": []}
    for i in range(len(test_Mags)):
        test_points["RAs"].extend(test_RAs)
        test_points["Decs"].extend(test_Decs)
        test_points["Mags"].extend([test_Mags[i]]*len(test_RAs))
        for j in range(len(test_RAs)):
            test_points["Observations"].append({})
            for filt in WFI_filt_fns:
                test_points["Observations"][-1][filt] = []

    for key in ["RAs", "Decs"]:
        test_points[key] = array(test_points[key])


    SN_data = {"nsne": nsne, "SN_table": Table([redshifts, RAs, Decs, daymaxes, [survey_fields]*nsne],
                                               names = ["redshifts", "RAs", "Decs", "daymaxes", "survey_fields"],
                                               dtype= ("f8", "f8", "f8", "f8", "S20")), "SN_observations": [],
               "CC_table": Table([redshifts, RAs, Decs, daymaxes, [survey_fields]*nsne],
                                 names = ["redshifts", "RAs", "Decs", "daymaxes", "survey_fields"],
                                 dtype= ("f8", "f8", "f8", "f8", "S20")), "SN_observations": [],
               "test_points": [test_points]}

    print("Getting SNCosmo models...")
    # I'm going to put the SN LC information into the per-SN dictionary list, as it needs to contain items like an SNCosmo model

    gal_backgrounds = make_galaxy_spectrum(redshifts)
    source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")

    for i in range(nsne):
        if i % 2000 == 0:
            print(i, nsne)
        MV, x1, c, host_mass, sncosmo_model = realize_SN(redshifts[i], daymaxes[i], source = source)
        

        SN_data["SN_observations"].append(dict(
            MV = MV, x1 = x1, c = c, host_mass = host_mass, sncosmo_model = sncosmo_model, gal_background = gal_backgrounds[i],
            dates = array([], dtype=float64), fluxes = array([], dtype=float64), dfluxes = array([], dtype=float64),
            filts = array([], dtype=(str, 10)), ispar = array([], dtype=bool), IFS_dates = [], IFS_fluxes = [], IFS_dfluxes = [], IFS_exptimes = [],
            found_date = None
        ))

    """
    dust = sncosmo.CCM89Dust()

    for i in range(1000):
        sncosmo_model = sncosmo.Model(source="nugent-sn2p", effects=[dust],
                                      effect_names=['host'],
                                      effect_frames=['rest'])

        SN_data["CC_observations"].append(dict(
            sncosmo_model = sncosmo_model, gal_background = gal_backgrounds[i],
            dates = array([], dtype=float64), fluxes = array([], dtype=float64), dfluxes = array([], dtype=float64),
            filts = array([], dtype=(str, 10)), ispar = array([], dtype=bool), IFS_dates = [], IFS_fluxes = [], IFS_dfluxes = [], IFS_exptimes = [],
            found_date = None
        ))
    """


    if salt2_model:
        if verbose:
            print("Observing ", tier_filters, tier_exptimes)


        SN_data = run_survey(SN_data = SN_data, square_degrees = square_degrees, tier_filters = tier_filters, tier_exptimes = tier_exptimes,
                             ground_depths = ground_depths, dithers = dithers, cadence = cadence,
                             fraction_time = fraction_time, total_survey_time = total_survey_time, hours_per_visit = hours_per_visit,
                             IFS_trigger_params = IFS_trigger_params, parallel_filters = parallel_filters, survey_duration = survey_duration)
    return SN_data

def merge_SN_data(SN_data, this_SN_data):

    if SN_data == {}:
        return this_SN_data
    else:

        assert list(SN_data.keys()) == list(this_SN_data.keys()), str(list(SN_data.keys())) + "_" + str(list(this_SN_data.keys()))

        print("Merging total_time_used")
        SN_data["total_time_used"] += this_SN_data["total_time_used"]

        print("Merging nsne")
        SN_data["nsne"] += this_SN_data["nsne"]

        print("Merging SN_observations")
        SN_data["SN_observations"].extend(this_SN_data["SN_observations"])

        print("Merging time_remaining_values")
        SN_data["time_remaining_values"].extend(this_SN_data["time_remaining_values"])

        print("Merging SN_table")
        SN_data["SN_table"] = vstack([SN_data["SN_table"], this_SN_data["SN_table"]])

        print("Merging observation_table")
        SN_data["observation_table"] = vstack([SN_data["observation_table"], this_SN_data["observation_table"]])

        print("Merging ordered_sol")
        SN_data["ordered_sol"] = vstack([SN_data["ordered_sol"], this_SN_data["ordered_sol"]])

        SN_data["test_points"].extend(this_SN_data["test_points"])

    return SN_data

################################################### Starting here ###################################################
master_zp = 25.0 # AB mag
picklefl = sys.argv[2]


survey_parameters = read_csv(sys.argv[1])

print("Reading PSFs...")
PSFs = initialize_PSFs(pixel_scales = [10, 15, 22], slice_scales = [30, 30, 22], PSF_source = survey_parameters["PSFs"])
PSFs_WFC = initialize_PSFs(pixel_scales = [10, 15, 22] + [60]*(survey_parameters["WFI_PSFs"] == "Euclid_PSF"),
                           slice_scales = [30, 30, 22] + [60]*(survey_parameters["WFI_PSFs"] == "Euclid_PSF"), PSF_source = survey_parameters["WFI_PSFs"])


slew_fn = file_to_fn(wfirst_data_path + "/pixel-level/input/" + survey_parameters["slew_table"], col = 2)


    
WFI_args = {"PSFs": PSFs_WFC, "source_dir": wfirst_data_path + "/pixel-level/input",
            "pixel_scale": survey_parameters["WFI_pixel_scale"],
            "dark_current": survey_parameters["WFI_dark_current"],
            "read_noise_floor": survey_parameters["WFI_read_noise_floor"],
            "read_noise_white": survey_parameters["WFI_read_noise_white"],
            "IPC": survey_parameters["interpixel_capacitance"],
            "TTel": survey_parameters["telescope_temperature"],
            "zodi_fl": file_to_fn(wfirst_data_path + "/pixel-level/input/" + survey_parameters["zodiacal_background"]),
            "bad_pixel_rate": survey_parameters["bad_pixel_rate"],
            "waves": arange(6000., 22500.1, 25.)}

IFS_args = {"PSFs": PSFs, "source_dir": wfirst_data_path + "/pixel-level/input",
            "pixel_scale": survey_parameters["IFU_pixel_scale"],
            "slice_scale": survey_parameters["IFU_slice_in_pixels"]*survey_parameters["IFU_pixel_scale"],
            "dark_current": survey_parameters["IFU_dark_current"],
            "read_noise_floor": survey_parameters["IFU_read_noise_floor"],
            "white_noise": survey_parameters["IFU_read_noise_white"],
            "IFURfl": file_to_fn(wfirst_data_path + "/pixel-level/input/" + survey_parameters["IFU_resolution"]),
            "zodifl": file_to_fn(wfirst_data_path + "/pixel-level/input/" + survey_parameters["zodiacal_background"]),
            "effareafl": file_to_fn(wfirst_data_path + "/pixel-level/input/" + survey_parameters["IFU_effective_area"]),
            "min_wave": survey_parameters["IFU_min_wave"], "max_wave": survey_parameters["IFU_max_wave"],
            "offset_par": int(around(survey_parameters["IFU_pixel_scale"]*0.25/0.005)), "offset_perp": 0,
            "IPC": survey_parameters["interpixel_capacitance"], "TTel": survey_parameters["telescope_temperature"]}


IFS_args["waves"] = get_spec_with_err(exp_time = 100, redshift = 1., phase = 0, show_plots = 0, **IFS_args)["obs_waves"]

ground_filt_fns, ground_obslambs, ground_five_sigma_one_hour, WFI_filt_fns = init_ground(survey_parameters["grizY_30s_ground_depths"])

ground_args = {"waves": ground_obslambs, "filt_fns": ground_filt_fns}

################################################### Run! ###################################################


assert isclose(sum(survey_parameters["tier_parameters"]["tier_fraction_time"]), 1), sum(survey_parameters["tier_parameters"]["tier_fraction_time"])
assert survey_parameters["maximum_trigger_fraction"] <= 1

rates_fn_full = file_to_fn(wfirst_data_path + "/pixel-level/input/" + survey_parameters["SN_rates"]) # per 1e-4 year Mpc^3 (h = 0.7)
rates_fn = lambda x: rates_fn_full(x)*survey_parameters["maximum_trigger_fraction"]

IFS_trigger_params = dict(adjust_each_SN_exp_time = survey_parameters["adjust_each_SN_exp_time"],
                          normalization_wavelength_range = survey_parameters["normalization_wavelength_range"],
                          SNR_set = [survey_parameters["shallow_SNR"], survey_parameters["medium_SNR"], survey_parameters["deep_SNR"], survey_parameters["reference_SNR"]/sqrt(survey_parameters["number_of_reference_dithers"])],
                          phases = [item[1] for item in survey_parameters["spectra_depth_and_phase"]])
IFS_trigger_params["SNRs"] = [IFS_trigger_params["SNR_set"][["shallow", "medium", "deep"].index(item[0])] for item in survey_parameters["spectra_depth_and_phase"]]
IFS_trigger_params["number_of_reference_dithers"] = survey_parameters["number_of_reference_dithers"]

print("IFS_trigger_params", IFS_trigger_params)

redshift_step = 0.05
redshift_set = arange(0.125, 2.475 + redshift_step/10., redshift_step)

survey_parameters["uniform_IFS_exptimes"] = {}
source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")

for redshift in redshift_set:
    print("Getting median...")
    mBs, x1s, colors, masses = make_SALT2_params(size = 10000)

    if max(IFS_args["waves"]) > IFS_trigger_params["normalization_wavelength_range"][1]*(1 + redshift):
        sncosmo_model = get_SNCosmo_model(redshift = redshift, x1 = median(x1s), c = median(colors), MV = median(mBs - colors), daymax = 0, source = source)
        survey_parameters["uniform_IFS_exptimes"][redshift] = get_IFS_exptimes(redshift = redshift, sncosmo_model = sncosmo_model, daymax = 0,
                                                                               IFS_trigger_params = IFS_trigger_params, gal_flamb = make_galaxy_spectrum([redshift], median_not_random = True)[0])

IFS_trigger_params["uniform_IFS_exptimes"] = survey_parameters["uniform_IFS_exptimes"]

SN_data = {}


for i in range(len(survey_parameters["tier_parameters"]["tier_name"])):
    IFS_trigger_params["trigger_redshifts"] = survey_parameters["tier_parameters"]["trigger_redshifts"][i]
    IFS_trigger_params["trigger_fraction"] = survey_parameters["tier_parameters"]["trigger_fraction"][i]


    this_SN_data = make_SNe(square_degrees = survey_parameters["tier_parameters"]["square_degrees"][i],
                            cadence = survey_parameters["tier_parameters"]["cadence"][i],
                            hours_per_visit = survey_parameters["hours_per_visit"],
                            survey_duration = survey_parameters["survey_duration"], rates_fn = rates_fn,
                            redshift_set = redshift_set,
                            IFS_trigger_params = IFS_trigger_params,
                            tier_filters = survey_parameters["tier_parameters"]["filters"][i],
                            tier_exptimes = quantize_time(survey_parameters["tier_parameters"]["exp_times_per_dither"][i]),
                            dithers = survey_parameters["tier_parameters"]["dithers_per_filter"][i],
                            ground_depths = get_ground_depths(survey_parameters, tier = i),
                            parallel_filters = survey_parameters["tier_parameters"]["parallel_filters"][i],
                            fraction_time = survey_parameters["tier_parameters"]["tier_fraction_time"][i],
                            total_survey_time = survey_parameters["total_survey_time"],
                            redshift_step = 0.05,
                            salt2_model = True, verbose = True, phase_buffer = 20,
                            survey_fields = survey_parameters["tier_parameters"]["tier_name"][i])
    this_SN_data["SN_table"]["RAs"] += i*20. # Make sure tiers don't overlap
    this_SN_data["observation_table"]["RA"] += i*20. # Make sure tiers don't overlap

    for j in range(len(this_SN_data["SN_observations"])):
        this_SN_data["SN_observations"][j]["sncosmo_model"] = None # Can't pickle SNCosmo model!
        this_SN_data["SN_observations"][j]["gal_background"] = this_SN_data["SN_observations"][j]["gal_background"](IFS_args["waves"])  # Can't pickle SNCosmo model!

    SN_data = merge_SN_data(SN_data, this_SN_data)


SN_data["survey_parameters"] = survey_parameters
SN_data["IFC_waves"] = IFS_args["waves"]

pickle.dump(SN_data, open(picklefl, 'wb'))

print("Done!")
