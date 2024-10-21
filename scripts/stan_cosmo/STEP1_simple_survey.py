import numpy as np
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
sys.path.append(wfirst_path + "/scripts/stan_cosmo/")
sys.path.append(wfirst_path + "/scripts/pixel-level/")
sys.path.append(wfirst_path + "/scripts/host/")
from host_generator import make_galaxy_spectrum
from generator import make_SALT2_params
from pixel_level_ETC2 import initialize_PSFs, get_imaging_SN, solve_for_exptime, get_spec_with_err
#from FileRead import readcol
import subprocess
import sncosmo
import pickle as pickle
import multiprocessing as mp
from astropy.table import Table, vstack
import matplotlib.path as MPLpath
import copy as cp
import tempfile
import time
import tqdm
from scipy.stats import percentileofscore
import matplotlib.patches as patches
whoami = subprocess.getoutput("whoami")


def file_to_fn(fl, col = 1):
    vals = loadtxt(fl)
    x = vals[:,0]
    y = vals[:,col]

    if np.max(y) - np.min(y) == 0:
        return lambda x: y[0]
    else:
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
    for key in ["total_survey_time", "survey_duration", "hours_per_visit", "SN_rates",
                "total_FoV", "active_FoV", "sn_model", "SN_number_poisson",
                "slew_table", "zodiacal_background", "telescope_temperature", "PSFs", "WFI_PSFs", "interpixel_capacitance",
                "WFI_dark_current", "WFI_read_noise_floor", "WFI_read_noise_white", "WFI_pixel_scale",
                "IFU_min_wave", "IFU_max_wave", "IFU_effective_area", "IFU_resolution", "IFU_pixel_scale", "IFU_slice_in_pixels", "IFU_dark_current", "IFU_read_noise_floor", "IFU_read_noise_white", "bad_pixel_rate"]:
        survey_parameters[key] = read_from_lines(lines, key)

    for key in ["grizY_30s_ground_depths"]:
        survey_parameters[key] = read_from_lines(lines, key, islist = True)


    # Now read each tier
    n_tiers = 0
    while read_from_lines(lines, "tier_name", item_number = n_tiers + 1) != None:
        n_tiers += 1
    print("n_tiers ", n_tiers)
    

    single_keys = ["tier_name", "square_degrees", "max_z", "max_SNe", "tier_fraction_time"]
    list_keys = ["filters", "exp_times_per_dither", "dithers_per_filter", "cadences", "cadence_offsets"]
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
    AB_mag = -2.5*log10(sum(ground_filt_fns[filt](ground_obslambs)*np.abs(flamb_SN)*ground_obslambs)/
                        sum(ground_filt_fns[filt](ground_obslambs)*(0.108848062485/ground_obslambs**2.)*ground_obslambs))
    
    AB_mag2 = -2.5*log10(sum(ground_filt_fns[filt](ground_obslambs + 10.)*np.abs(flamb_SN)*ground_obslambs)/
                        sum(ground_filt_fns[filt](ground_obslambs + 10.)*(0.108848062485/ground_obslambs**2.)*ground_obslambs))
    return AB_mag, AB_mag2 - AB_mag


    

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
    for filt in ["R062", "Z087", "Y106", "J129", "H158", "F184", "K213", "W146", "P100", "Euclid_Y", "Euclid_J", "Euclid_H"]:
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
            tfs = SN_data["SN_observations"][i]["true_fluxes"][inds]
            

            new_f = sum(ws * fs)/sum(ws)
            new_tf = sum(ws * tfs)/sum(ws)
            new_df = 1./sqrt(sum(ws))

            SN_data["SN_observations"][i]["fluxes"][inds[0]] = new_f
            SN_data["SN_observations"][i]["true_fluxes"][inds[0]] = new_tf
            SN_data["SN_observations"][i]["dfluxes"][inds[0]] = new_df
            SN_data["SN_observations"][i]["eminus_per_s"][inds[0]] = sum(SN_data["SN_observations"][i]["eminus_per_s"][inds] * ws)/sum(ws)
            SN_data["SN_observations"][i]["dmag_d10A"][inds[0]] = sum(SN_data["SN_observations"][i]["dmag_d10A"][inds] * ws)/sum(ws)
            SN_data["SN_observations"][i]["ispar"][inds[0]] = any(SN_data["SN_observations"][i]["ispar"][inds])

            for key in ["dates", "true_fluxes", "fluxes", "dfluxes", "filts", "eminus_per_s", "dmag_d10A"]:
                SN_data["SN_observations"][i][key] = SN_data["SN_observations"][i][key][:nobs+1-len(inds)]

    return SN_data


def get_dmag_d10A(flamb_SN, obslambs, filt_fn):
    AB_mag = -2.5*log10(sum(filt_fn(obslambs)*np.abs(flamb_SN)*obslambs)/
                        sum(filt_fn(obslambs)*(0.108848062485/obslambs**2.)*obslambs))
    
    AB_mag2 = -2.5*log10(sum(filt_fn(obslambs + 10.)*np.abs(flamb_SN)*obslambs)/
                         sum(filt_fn(obslambs + 10.)*(0.108848062485/obslambs**2.)*obslambs))

    return AB_mag2 - AB_mag
    


def imaging_ETC_wrapper(SN_data, ind, row_to_add, current_date):
    #print("row_to_add", row_to_add)
    #print("current_date", current_date)
    


    
    f_lamb_SN = SN_data["SN_observations"][ind]["sncosmo_model"].flux(current_date, WFI_args["waves"])
    WFI_args["mdl"] = f_lamb_SN
                
    ETC_result = get_imaging_SN(redshift = 0, exp_time = row_to_add["exptime"], gal_flamb = SN_data["SN_observations"][ind]["gal_background"],
                                effective_meters2_fl = WFI_filt_fns[row_to_add["filt"]], phase = 0,
                                offset_par = random.random()*22, offset_perp = random.random()*22, **WFI_args)
    
    #ETC_result = {"AB_mag": 24, "PSF_phot_S/N": 5.}
    #print("ETC_result", ETC_result.keys())
    # ETC_result dict_keys(['PSF_phot_S/N', 'SN_photons/s', 'PSF-weighted total/s', 'AB_mag', 'zodi/s', 'thermal/s'])

    
    AB_flux = 10.**(-0.4*(ETC_result["AB_mag"] - master_zp))
    SNR_total = ETC_result["PSF_phot_S/N"]
    
    flux_with_noise = AB_flux + random.normal()*AB_flux/SNR_total
    SN_data["SN_observations"][ind]["dates"] = append(SN_data["SN_observations"][ind]["dates"], current_date)
    SN_data["SN_observations"][ind]["true_fluxes"] = append(SN_data["SN_observations"][ind]["true_fluxes"], AB_flux)
    # eminus_per_s = array([], dtype=float64), dmag_d10A
    SN_data["SN_observations"][ind]["eminus_per_s"] = append(SN_data["SN_observations"][ind]["eminus_per_s"], ETC_result["PSF-weighted total/s"])
    SN_data["SN_observations"][ind]["dmag_d10A"] = append(SN_data["SN_observations"][ind]["dmag_d10A"], get_dmag_d10A(flamb_SN = f_lamb_SN,
                                                                                                                      obslambs = WFI_args["waves"],
                                                                                                                      filt_fn = WFI_filt_fns[row_to_add["filt"]]))
    SN_data["SN_observations"][ind]["fluxes"] = append(SN_data["SN_observations"][ind]["fluxes"], flux_with_noise)
    SN_data["SN_observations"][ind]["dfluxes"] = append(SN_data["SN_observations"][ind]["dfluxes"], AB_flux/SNR_total)
    SN_data["SN_observations"][ind]["filts"] = append(SN_data["SN_observations"][ind]["filts"], row_to_add["filt"])
    SN_data["SN_observations"][ind]["ispar"] = append(SN_data["SN_observations"][ind]["ispar"], row_to_add["SNind"] > -1)

    return SN_data


def prism_ETC_wrapper(SN_data, ind, row_to_add, current_date):
    f_lamb_SN = SN_data["SN_observations"][ind]["sncosmo_model"].flux(current_date, IFS_args["waves"])

    assert all(1 - isnan(f_lamb_SN))
    assert row_to_add["filt"] == "P100"
    
    #args = cp.deepcopy(IFS_args)
    IFS_args["mdl"] = f_lamb_SN
    
    ETC_result = get_spec_with_err(redshift = 0, exp_time = row_to_add["exptime"], gal_flamb = lambda x:0, bad_pixel_rate = survey_parameters["bad_pixel_rate"],
                                   show_plots = 0, **IFS_args)
    
    #print("ETC_result", ETC_result.keys(), ETC_result["PSF_wghtd_e_allbutdark"], "PSF_wghtd_e_allbutdark")
    #fldkajfldskj
    SN_data["SN_observations"][ind]["IFS_dates"].append(current_date)
    SN_data["SN_observations"][ind]["IFS_exptimes"].append(row_to_add["exptime"])


    SN_data["SN_observations"][ind]["IFS_eminus_per_s"].append(ETC_result["PSF_wghtd_e_allbutdark"]/row_to_add["exptime"] + 1.2)
    
    noise = ETC_result["f_lamb_SN"]/ETC_result["spec_S/N"]
    assert all(1 - isnan(noise) - isinf(noise))
    SN_data["SN_observations"][ind]["IFS_true_fluxes"].append(ETC_result["f_lamb_SN"])
    
    SN_data["SN_observations"][ind]["IFS_fluxes"].append(
        ETC_result["f_lamb_SN"] + random.normal(size = len(ETC_result["f_lamb_SN"]))*noise)
    SN_data["SN_observations"][ind]["IFS_dfluxes"].append(noise)
    return SN_data


def run_observation_through_ETC(SN_data, row_to_add, current_date):
    print("Observing", row_to_add, current_date)
    
    phases = (current_date - SN_data["SN_table"]["daymaxes"])/(1. + SN_data["SN_table"]["redshifts"])

    active_SN_mask = (phases >= -15.)*(phases <= 45.)
    if any(active_SN_mask):
        inds = where(active_SN_mask)[0]
        

        for ind in inds:
            #args = cp.deepcopy(WFI_args)
            if row_to_add["filt"] != "P100":
                if np.random.random() < survey_parameters["active_FoV"]/survey_parameters["total_FoV"]: # Approximate fill factor:
                    SN_data = imaging_ETC_wrapper(SN_data = SN_data, ind = ind, row_to_add = row_to_add, current_date = current_date)
            else:
                if np.random.random() <	survey_parameters["active_FoV"]/survey_parameters["total_FoV"]:
                    SN_data = prism_ETC_wrapper(SN_data = SN_data, ind = ind, row_to_add = row_to_add, current_date = current_date)


    for ind in range(len(SN_data["test_points"][0]["Mags"])):
        f_lamb_SN = 0.10884806248*10.**(-0.4*SN_data["test_points"][0]["Mags"][ind]) /(WFI_args["waves"]**2.)

        WFI_args["mdl"] = f_lamb_SN

        ETC_result = get_imaging_SN(redshift = 0, exp_time = row_to_add["exptime"], gal_flamb = lambda x:0,
                                    effective_meters2_fl = WFI_filt_fns[row_to_add["filt"]], phase = 0,
                                    offset_par = random.random()*22, offset_perp = random.random()*22, **WFI_args)

        SNR_total = ETC_result["PSF_phot_S/N"]
        SN_data["test_points"][0]["Observations"][ind][row_to_add["filt"]].append(SNR_total)
    
    return SN_data


def run_observation_through_ground_ETC(SN_data, row_to_add, current_date, ground_depths):
    phases = (current_date - SN_data["SN_table"]["daymaxes"])/(1. + SN_data["SN_table"]["redshifts"])

    active_SN_mask = (phases >= -15.)*(phases <= 45.)
    
    inds = where(active_SN_mask)[0]
    
    for ind in inds:
        f_lamb_SN = SN_data["SN_observations"][ind]["sncosmo_model"].flux(current_date, ground_args["waves"])
        
        AB_mag, dAB_mag_d10A = get_ground_AB_mag(f_lamb_SN, row_to_add["filt"])        
        
        AB_flux = 10.**(-0.4*(AB_mag - master_zp))
        AB_err = 0.2*10.**(-0.4*(ground_depths[row_to_add["filt"]] - master_zp))
        SNR_total = AB_flux/AB_err

    
        flux_with_noise = AB_flux + random.normal()*AB_flux/SNR_total
        SN_data["SN_observations"][ind]["dates"] = append(SN_data["SN_observations"][ind]["dates"], current_date)
        SN_data["SN_observations"][ind]["true_fluxes"] = append(SN_data["SN_observations"][ind]["true_fluxes"], flux_with_noise)
        SN_data["SN_observations"][ind]["eminus_per_s"] = append(SN_data["SN_observations"][ind]["eminus_per_s"], 10000.) # Ground-based data has no CRNL
        SN_data["SN_observations"][ind]["dmag_d10A"] = append(SN_data["SN_observations"][ind]["dmag_d10A"], dAB_mag_d10A)
        SN_data["SN_observations"][ind]["fluxes"] = append(SN_data["SN_observations"][ind]["fluxes"], flux_with_noise)
        SN_data["SN_observations"][ind]["dfluxes"] = append(SN_data["SN_observations"][ind]["dfluxes"], AB_flux/SNR_total)
        SN_data["SN_observations"][ind]["filts"] = append(SN_data["SN_observations"][ind]["filts"], row_to_add["filt"])


        
    return SN_data


def quantize_time(t):
    return np.ceil(np.array(t)/3.04001)*3.04 # The ...01 is for numerical accuracy.

def get_slew_time(square_degrees):
    pointings = square_degrees/survey_parameters["total_FoV"]
    return pointings*55.

def add_observations_to_sndata(SN_data, rows_to_add, current_date, ground_depths, square_degrees):
    rows_added = Table(names = ["date", "filt", "exptime", "RA", "dec", "orient", "instr", "SNind", "tier"],
                       dtype= ("f8", "S10", "f8", "f8", "f8", "f8", "S10", "i4", "S40"))

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
            tmp_slew_time = get_slew_time(square_degrees)

            
            time_used += tmp_slew_time
            slew_time += tmp_slew_time
        
        for ind in inds:
            # For this filter, for current_date, for row "ind"
            if rows_to_add["instr"][ind] != "ground":
                time_used += quantize_time(rows_to_add["exptime"][ind])*(square_degrees/survey_parameters["total_FoV"])
                SN_data = run_observation_through_ETC(SN_data, rows_to_add[ind], current_date)
            else:
                #print("ground_depths", ground_depths)
                SN_data = run_observation_through_ground_ETC(SN_data, rows_to_add[ind], current_date, ground_depths)

            rows_added.add_row(rows_to_add[ind])
            rows_to_add.remove_row(ind) # Done with this observation!

        # Finished adding all data for this filter, now, take variance-weighted mean of all repeat observations for a given SN (either because of dithers, or parallels)
        SN_data = compress_lc_data(SN_data, current_date, filt)

    return SN_data, rows_to_add, rows_added, time_used, slew_time


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


def plan_and_add_cadence(SN_data, wfi_time_left, next_date, rows_to_add,
                         tier_cadences, tier_cadence_offsets,
                         tier_filters, tier_exptimes, tier_dithers, tier):


    for filt, expt, dith, cadence, cadence_offset in zip(tier_filters, tier_exptimes, tier_dithers, tier_cadences, tier_cadence_offsets):
        dates_for_this_filt = np.arange(cadence_offset, 10000, cadence)
        if len(filt) > 1:
            # If WFI, not ground
            if wfi_time_left:
                for j in range(dith):
                    if np.min(np.abs(dates_for_this_filt - next_date)) < 0.001:
                        rows_to_add.add_row((next_date, filt, expt,
                                             0.0, # RA
                                             0.0, # Dec
                                             0.0, # Angle
                                             "WFI", -1, tier))


            else:
                pass # Out of time, nothing to add
        else:
            # Cadence offsets don't work here
            assert cadence_offset == 0
            rows_to_add.add_row((next_date, filt, 0, 0, 0, 0, "ground", -1, tier))

    return rows_to_add


def run_survey(SN_data, square_degrees, tier_filters, tier_exptimes, tier, ground_depths, dithers, tier_cadences, tier_cadence_offsets, total_survey_time, hours_per_visit, survey_duration):

    starting_time = 31557000*total_survey_time  # Seconds in total_survey_time years
    total_time_left = starting_time
    total_slew_time = 0.0

    # Start with empty table
    # instr is WFI or ground
    # SNind is index of SN (if any) in IFS
    # RA, dec are WFI, not IFS

    observation_table = Table(names = ["date", "filt", "exptime", "RA", "dec", "orient", "instr", "SNind", "tier"],
                              dtype= ("f8", "S10", "f8", "f8", "f8", "f8", "S10", "i4", "S40"))

    # RA is x, Dec is y

    print(observation_table)


    # Start with empty table
    rows_to_add = Table(names = ["date", "filt", "exptime", "RA", "dec", "orient", "instr", "SNind", "tier"],
                        dtype= ("f8", "S10", "f8", "f8", "f8", "f8", "S10", "i4", "S40"))



    possible_cadence_steps = []
    for i in range(len(tier_cadences)):
        possible_cadence_steps += list(np.arange(tier_cadence_offsets[i], 365.24*survey_duration + 100, tier_cadences[i]))
    possible_cadence_steps = list(set(possible_cadence_steps))
    possible_cadence_steps.sort()

    print(possible_cadence_steps)
    
    
    assert possible_cadence_steps[0] == 0
    current_date = possible_cadence_steps[0]
    SN_data["time_remaining_values"] = [[(current_date, total_time_left, total_time_left)]]

    while (len(rows_to_add["filt"]) > 0 or current_date == 0) and (current_date < 365.24*survey_duration) and (total_time_left >= 0): # The >= 0 is important for ground-based tiers with no Roman

        # Step 1: Add previously planned observations
        SN_data, rows_to_add, rows_added, time_used, slew_time = add_observations_to_sndata(SN_data = SN_data, rows_to_add = rows_to_add,
                                                                                            current_date = current_date, square_degrees = square_degrees,
                                                                                            ground_depths = ground_depths)
        for i in range(len(rows_added)):
            observation_table.add_row(rows_added[i])
        total_time_left -= time_used
        total_slew_time += slew_time
        
        # Step 2: Find SNe in observations
        SN_data = find_SNe(SN_data, current_date)
        

        # Step 3: Plan cadenced observations for next time

        estimated_time_left = total_time_left - sum(rows_to_add["exptime"] + 50.)
        wfi_time_left = estimated_time_left > 0
        # This is an approximate ending condition; we'll see how close total_time_left comes to zero.
            
        rows_to_add = plan_and_add_cadence(SN_data = SN_data, wfi_time_left = wfi_time_left, next_date = possible_cadence_steps[1],
                                           rows_to_add = rows_to_add,
                                           tier_cadences = tier_cadences, tier_cadence_offsets = tier_cadence_offsets,
                                           tier_filters = tier_filters, tier_exptimes = tier_exptimes, tier_dithers = dithers, tier = tier)
            

        print(current_date, "rows_to_add", rows_to_add)
        
        print("\nEnd of day", current_date, "time left", total_time_left, time.asctime(), '\n')

        del possible_cadence_steps[0]
        current_date = possible_cadence_steps[0]
        SN_data["time_remaining_values"][0].append((current_date, total_time_left, estimated_time_left))
    
        print("date ", current_date, "total_time_left ", total_time_left)


    print("Done with tier. Slew time = ", total_slew_time)
    SN_data["total_time_used"] = starting_time - total_time_left
    SN_data["observation_table"] = observation_table

    return SN_data


def make_random_positions(inner_radius, outer_radius, nsne):
    CDF = random.random(size = nsne)
    rs = np.sqrt((1. - CDF)*inner_radius**2. + CDF*outer_radius**2.)
    thetas = random.random(size = nsne)*2*pi
        
    RAs = rs*cos(thetas)
    Decs = rs*sin(thetas)
    return RAs, Decs


def make_SNe(square_degrees, tier_cadences, tier_cadence_offsets, survey_duration, hours_per_visit,
             rates_fn, redshift_set, tier_filters, tier_exptimes, ground_depths, dithers, SN_number_poisson,
             total_survey_time, max_z, max_SNe, redshift_step = 0.05, salt2_model = True, verbose = False, phase_buffer = 20, survey_fields = "None"):
    #assert square_degrees <= 5000, "Should use more accurate formula for large surveys!"


    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    frac_of_sky = (square_degrees)/(4.*pi*(180./pi)**2.)
    
    if verbose:
        print("frac_of_sky ", frac_of_sky)
        

    volume_in_shells = cosmo.comoving_volume(redshift_step/2. + redshift_set).value - cosmo.comoving_volume(redshift_set - redshift_step/2.).value
    SNe_in_survey_field = volume_in_shells*frac_of_sky*rates_fn(redshift_set) * survey_duration/(1. + redshift_set) * 1e-4

    if verbose:
        print(SNe_in_survey_field, sum(SNe_in_survey_field))

    print("SN_number_poisson", SN_number_poisson)
    
    if SN_number_poisson:
        SNe_actual = [random.poisson(SNe_in_survey_field[i]) for i in range(len(redshift_set)) if redshift_set[i] <= max_z]
    else:
        SNe_actual = [int(np.around(SNe_in_survey_field[i])) for i in range(len(redshift_set)) if redshift_set[i] <= max_z]

    if verbose:
        print(SNe_actual, sum(SNe_actual))
    redshifts = []
    for i in range(len(redshift_set)):
        if redshift_set[i] <= max_z:
            redshifts += [redshift_set[i]]*SNe_actual[i]
    redshifts = array(redshifts)

    assert len(redshifts) == sum(SNe_actual)

    daymaxes = random.random(size = sum(SNe_actual))*survey_duration*365
    not_at_edge = where((daymaxes/(1. + redshifts) > phase_buffer)*((survey_duration*365 - daymaxes)/(1. + redshifts) > phase_buffer))

    redshifts = redshifts[not_at_edge]
    daymaxes = daymaxes[not_at_edge]

    if len(redshifts) > max_SNe:
        the_mask = np.array([1] * max_SNe + [0] * (len(redshifts) - max_SNe))
        np.random.shuffle(the_mask)

        inds = np.where(the_mask)
        redshifts = redshifts[inds]
        daymaxes = daymaxes[inds]
        
    
    nsne = len(redshifts)

    if verbose:
        print("SNe not near beginning/end: ", nsne)

    outer_radius = np.sqrt(square_degrees/np.pi)
    RAs, Decs = make_random_positions(inner_radius = 0.0, outer_radius = outer_radius, nsne = nsne)
    test_RAs, test_Decs = make_random_positions(inner_radius = 0.0, outer_radius = outer_radius, nsne = 10)
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
                                               dtype= ("f8", "f8", "f8", "f8", "S40")), "SN_observations": [],
               "CC_table": Table([redshifts, RAs, Decs, daymaxes, [survey_fields]*nsne],
                                 names = ["redshifts", "RAs", "Decs", "daymaxes", "survey_fields"],
                                 dtype= ("f8", "f8", "f8", "f8", "S40")), "SN_observations": [],
               "test_points": [test_points]}

    print("Getting SNCosmo models...")
    # I'm going to put the SN LC information into the per-SN dictionary list, as it needs to contain items like an SNCosmo model

    gal_backgrounds = make_galaxy_spectrum(redshifts)

    if survey_parameters["sn_model"] == "salt2_extended":
        source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")
    elif survey_parameters["sn_model"] == "SALT3.NIR_WAVEEXT":
        source = sncosmo.SALT3Source(modeldir=wfirst_data_path + "/SALT3.NIR_WAVEEXT/")
    else:
        assert 0, "Unknown sn_model " + survey_parameters["sn_model"]
        
    for i in tqdm.trange(nsne):
        MV, x1, c, host_mass, sncosmo_model = realize_SN(redshifts[i], daymaxes[i], source = source)

        if i == 0:
            print(sncosmo_model)

        SN_data["SN_observations"].append(dict(
            MV = MV, x1 = x1, c = c, host_mass = host_mass, sncosmo_model = sncosmo_model, gal_background = gal_backgrounds[i],
            dates = array([], dtype=float64), true_fluxes = array([], dtype=float64), fluxes = array([], dtype=float64), dfluxes = array([], dtype=float64),
            filts = array([], dtype=(str, 10)), ispar = array([], dtype=bool), IFS_dates = [], IFS_fluxes = [], IFS_true_fluxes = [], IFS_dfluxes = [], IFS_exptimes = [], IFS_eminus_per_s = [],
            found_date = None, eminus_per_s = array([], dtype=float64), dmag_d10A = array([], dtype=float64)
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


    SN_data = run_survey(SN_data = SN_data, square_degrees = square_degrees, tier_filters = tier_filters, tier_exptimes = tier_exptimes, tier = survey_fields,
                         ground_depths = ground_depths, dithers = dithers, tier_cadences = tier_cadences, tier_cadence_offsets = tier_cadence_offsets,
                         total_survey_time = total_survey_time, hours_per_visit = hours_per_visit, survey_duration = survey_duration)
    
            
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

        SN_data["test_points"].extend(this_SN_data["test_points"])

    return SN_data

################################################### Starting here ###################################################

if __name__ == "__main__":
    master_zp = 25.0 # AB mag
    picklefl = sys.argv[2]


    survey_parameters = read_csv(sys.argv[1])

    print("Reading PSFs...")
    PSFs = initialize_PSFs(pixel_scales = [10, 15, 22, 22], slice_scales = [30, 30, 22, 110], PSF_source = survey_parameters["PSFs"])
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




    rates_fn_full = file_to_fn(wfirst_data_path + "/pixel-level/input/" + survey_parameters["SN_rates"]) # per 1e-4 year Mpc^3 (h = 0.7)
    rates_fn = rates_fn_full


    redshift_step = 0.05
    max_max_z = max(survey_parameters["tier_parameters"]["max_z"])
    print("max_max_z", max_max_z)
    
    redshift_set = arange(0.075, max_max_z + redshift_step/10., redshift_step)

    #source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")

    SN_data = {}


    inner_radius = 0.

    for i in range(len(survey_parameters["tier_parameters"]["tier_name"])):
        print("Running tier", survey_parameters["tier_parameters"]["tier_name"][i])


        this_SN_data = make_SNe(square_degrees = survey_parameters["tier_parameters"]["square_degrees"][i],
                                tier_cadences = survey_parameters["tier_parameters"]["cadences"][i],
                                tier_cadence_offsets = survey_parameters["tier_parameters"]["cadence_offsets"][i],
                                hours_per_visit = survey_parameters["hours_per_visit"],
                                survey_duration = survey_parameters["survey_duration"], rates_fn = rates_fn,
                                redshift_set = redshift_set,
                                SN_number_poisson = survey_parameters["SN_number_poisson"],
                                tier_filters = survey_parameters["tier_parameters"]["filters"][i],
                                tier_exptimes = quantize_time(survey_parameters["tier_parameters"]["exp_times_per_dither"][i]),
                                dithers = survey_parameters["tier_parameters"]["dithers_per_filter"][i],
                                ground_depths = get_ground_depths(survey_parameters, tier = i),
                                total_survey_time = survey_parameters["total_survey_time"]*survey_parameters["tier_parameters"]["tier_fraction_time"][i],
                                max_z = survey_parameters["tier_parameters"]["max_z"][i],
                                max_SNe = survey_parameters["tier_parameters"]["max_SNe"][i],
                                redshift_step = 0.05,
                                salt2_model = True, verbose = True, phase_buffer = 20,
                                survey_fields = survey_parameters["tier_parameters"]["tier_name"][i])


        for j in range(len(this_SN_data["SN_observations"])):
            this_SN_data["SN_observations"][j]["sncosmo_model"] = None # Can't pickle SNCosmo model!
            print('this_SN_data["SN_observations"][j]["gal_background"]', this_SN_data["SN_observations"][j]["gal_background"])
            
            this_SN_data["SN_observations"][j]["gal_background"] = this_SN_data["SN_observations"][j]["gal_background"](IFS_args["waves"][::10])  # Can't pickle SNCosmo model!

        SN_data = merge_SN_data(SN_data, this_SN_data)


    SN_data["survey_parameters"] = survey_parameters
    SN_data["IFC_waves"] = IFS_args["waves"]

    pickle.dump(SN_data, open(picklefl, 'wb'))

    print("Done!")
