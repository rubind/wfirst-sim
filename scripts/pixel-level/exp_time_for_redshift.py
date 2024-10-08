from copy import deepcopy
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pixel_level_ETC2 import initialize_PSFs, get_imaging_SN, solve_for_exptime, get_spec_with_err, flamb_to_photons_per_wave, resolution_to_wavelengths
import sys
from FileRead import readcol, writecol, format_file
import subprocess
import argparse
import glob
import sncosmo
from astropy import cosmology
import sys
import os
wfirst_path = os.environ["WFIRST"]
wfirst_data_path = os.environ["WFIRST_SIM_DATA"]
sys.path.append(wfirst_path + "/scripts/host/")
from host_generator import make_galaxy_spectrum
sys.path.append(wfirst_path + "/scripts/stan_cosmo/")
from generator import make_SALT2_params
from scipy.stats import scoreatpercentile
from DavidsNM import miniNM_new
import argparse
import tqdm

def file_to_fn(fl, col = 1):
    vals = np.loadtxt(fl)
    x = vals[:,0]
    y = vals[:,col]

    return interp1d(x, y, kind = 'linear')

def get_zodi(eff_area_fl, pixel_scale = 0.11):
    dwaves = 1.
    waves = np.arange(4000., 25000. + dwaves, dwaves)

    log10zodi_fn = file_to_fn(wfirst_data_path + "/pixel-level/input/aldering.txt")
    m2fn = file_to_fn(eff_area_fl)
    zodi_fn = lambda x: 10.**(log10zodi_fn(x))

    return sum(flamb_to_photons_per_wave(zodi_fn(waves), m2fn(waves), waves, dwaves))*pixel_scale**2.



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
    x1 = np.random.normal() - 0.25
    c = np.random.normal()*0.1 + np.random.exponential()*0.1 - 0.1*np.log(2.)
    host_mass = np.random.normal() + 10.

    MV = -19.08 + np.random.normal()*0.12 + np.random.normal()*0.055*redshift - 0.13*x1 + 2.1*c - (host_mass > 10.)*0.08 # beta is 2.1, as this is MV
    """
    
    MB, x1, color, mass = make_SALT2_params(size = 1)
    [MB, x1, color, mass] = [item[0] for item in [MB, x1, color, mass]]

    MB += np.random.normal()*0.1 + 0.055*redshift*np.random.normal() + 5/np.log(10.)*(0.001/redshift)*np.random.normal()

    sncosmo_model = get_SNCosmo_model(redshift = redshift,
                                      x1 = x1,
                                      c = color,
                                      MV = MB - color,
                                      daymax = daymax, source = source)


    return MB - color, x1, color, mass, sncosmo_model

def run_ETC(redshift, phase, source, exp_time, gal_flambs):
    SNRs = []

    for i in range(nsne):
        NA, NA, NA, NA, sncosmo_model = realize_SN(redshift = redshift, daymax = 0, source = source)

        f_lamb_SN = sncosmo_model.flux(phase*(1. + redshift), WFI_args["waves"])

        
        ETC_result = get_imaging_SN(redshift = 0, exp_time = exp_time, gal_flamb = gal_flambs[i],
                                    effective_meters2_fl = effective_meters2_fl, phase = 0, mdl = f_lamb_SN,
                                    offset_par = int(np.around(np.random.random()*22)), offset_perp = int(np.around(np.random.random()*22)), **WFI_args)
        SNRs.append(ETC_result["PSF_phot_S/N"])

    return SNRs

def run_ETC_prism(redshift, source, exp_time):
    inds = np.where((WFI_args["waves"]/(1. + redshift) >= 5000)*(WFI_args["waves"]/(1. + redshift) <= 6000))
    SNRs = []

    for i in range(nsne):
        NA, NA, NA, NA, sncosmo_model = realize_SN(redshift = redshift, daymax = 0, source = source)

        SNR_per_epoch = []
        for phase in np.arange(-15, 45, 5./(1. + redshift)):
            f_lamb_SN = sncosmo_model.flux(phase*(1. + redshift), WFI_args["waves"])
            #plt.figure(5)
            #plt.plot(WFI_args["waves"], f_lamb_SN, label = str(phase))

            ETC_result = get_spec_with_err(exp_time = exp_time, mdl = f_lamb_SN,
                                           offset_par = int(np.around(np.random.random()*22)), offset_perp = int(np.around(np.random.random()*22)), **WFI_args)

            this_StoN = ETC_result["spec_S/N"]

            SNR_per_epoch.append(np.sqrt(sum(this_StoN[inds]**2.)))
        #plt.legend(loc = 'best')
        ##plt.savefig("tmp.pdf")
        #plt.close()
        
        SNRs.append(float(np.sqrt(np.dot(SNR_per_epoch, SNR_per_epoch))))

    return SNRs


def make_SN_curve(P, exp_time):
    readouts = exp_time/3.04
    return exp_time/np.sqrt(P[0] + P[1]/exp_time*(readouts - 1)/(readouts + 1) + P[2]*exp_time)

def chi2fn(P, thedata):
    exp_times, SNs = thedata[0]
    
    if any(P < 0):
        return 1e100

    model = make_SN_curve(P, exp_times)
    resid = (SNs - model)/SNs
    return np.dot(resid, resid)

def fit_curve(exp_times, SNs):
    bestF = 1e100

    for ministart in ([1., 1., 1.],
                      [10., 10., 10.],
                      [0.2, 0.2, 0.2],
                      [50., 50., 50.]):
        

        P, F, NA = miniNM_new(ministart = ministart, miniscale = [0.1, 0.1, 0.1], chi2fn = chi2fn, passdata = [exp_times, SNs], compute_Cmat = False)
    
        if F < bestF:
            bestF = F
            bestP = P

    return make_SN_curve(bestP, exp_times)

def SNR_label(exp_times, curve, SNRtargets):
    label = ""
    vals = {}

    for SNRtarget in SNRtargets:
        ifn = interp1d(curve, exp_times, kind = 'linear', bounds_error = False, fill_value = 0)
        label += "S/N %.1f: %.2g s " % (SNRtarget, ifn(SNRtarget))
        vals[SNRtarget] = ifn(SNRtarget)
    return label, vals


parser = argparse.ArgumentParser()

parser.add_argument("--filt", help="Filter to Run", type=str)
parser.add_argument("--ttel", help="Telescope Temperature", type=float, default = 264)
parser.add_argument("--wfcdark", help="WFC dark current", type=float, default = 0.015)
parser.add_argument("--wfcrnf", help="WFC read-noise floor", type=float, default = 5.)
parser.add_argument("--wfcrnw", help="WFC white read noise", type=float, default = 20.)
parser.add_argument("--zodi", help="Zodi File", type=str, default = "aldering")
parser.add_argument("--zlots", help="Lots of Redshifts", type=float, default = 0.)
parser.add_argument("--nsne", help="Number of SNe", type=int, default = 400)
parser.add_argument("--SNRpeaktarg", help="S/N per point at peak", type=float, default = 10)
parser.add_argument("--SNRpeakred", help="S/N per point at this redshift", type=float, default = 1)
parser.add_argument("--addhost", type=int, default = 1)
parser.add_argument("--SNmodel", type=str)


opts = parser.parse_args()

print(opts)

if opts.filt == "P100":
    prism = 1
else:
    prism = 0

print("Reading PSFs...")

if prism:
    psfnm = "Prism_3Element_PSF_SCA08"
    PSFs = initialize_PSFs(pixel_scales = [22], slice_scales = [100], PSF_source = psfnm)

    zodi_per_pix = get_zodi(eff_area_fl = wfirst_data_path + "/pixel-level/input/P100.txt")
    zodi_per_pix *= 1.05 # For host-galaxy light
    print("zodi_per_pix ", zodi_per_pix)

    waves, dwaves = resolution_to_wavelengths(source_dir = wfirst_data_path, IFURfl = file_to_fn(wfirst_data_path + "/pixel-level/input/Prism_R.txt"), min_wave = 7500., max_wave = 18000.)


    WFI_args = {"gal_flamb": lambda x:0, "pixel_scale": 0.11, "redshift": 0, "phase": 0, "white_noise": 30, "read_noise_floor": 5., 
                "slice_scale": 0.5, "dark_current": zodi_per_pix + 0.015, "PSFs": PSFs, "waves": waves,
                "effareafl": file_to_fn(wfirst_data_path + "/pixel-level/input/P100.txt"), "TTel": 200,
                "restframe_bins": (5000, 6000), "psfsize": 5, "zodifl": lambda x: -20}

else:
    psfnm = "WebbPSF_WFC"
    
    PSFs = initialize_PSFs(pixel_scales = [22], slice_scales = [22], PSF_source = psfnm)
    WFI_args = {"PSFs": PSFs, "source_dir": wfirst_data_path + "/pixel-level/input/",
                "pixel_scale": 0.11, "dark_current": opts.wfcdark, "read_noise_white": opts.wfcrnw, "read_noise_floor": opts.wfcrnf,
                "IPC": 0.02,
                "TTel": opts.ttel,
                "zodi_fl": file_to_fn(wfirst_data_path + "/pixel-level/input/" + opts.zodi + ".txt"),
                "bad_pixel_rate": 0.01,
                "waves": np.arange(4500., 25000.1, 25.)}
    

nsne = opts.nsne

if 1:
    if prism:
        if opts.zlots == 0:
            redshifts = np.arange(0.3, 2.01, 0.1)
        else:
            redshifts = [float(opts.zlots)]
        phases = [-99]
    else:
        if opts.zlots:
            redshifts = np.arange(0.2, 3.51, 0.1)
            phases = [0]
        else:
            redshifts = [0.4, 0.8, 1.0, 1.2, 1.7, 2.0]
            phases = [-10, -8, -6, 0]

    exp_times = 10**np.arange(0.5, 5.01, 0.05)
    if opts.SNmodel == "salt2_extended":
        source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")
    elif opts.SNmodel == "SALT3.NIR_WAVEEXT":
        source = sncosmo.SALT3Source(modeldir=wfirst_data_path + "/SALT3.NIR_WAVEEXT/")
    else:
        assert 0
        
    effective_meters2_fl = file_to_fn(wfirst_data_path + "/pixel-level/input/" + opts.filt + ".txt")

    for sqrtt in [0]:
        plt.figure(figsize = (6*len(redshifts), 4*len(phases)))
        pltname = "StoN_vs_exptime_%s_TTel_%.1f_%s%s_%s_dark=%.3f_rnf=%.1f_rnw=%.1f_SNRtarg=%.2f@%.2f_addhost=%i_mod=%s_zlots=%.1f" % (opts.filt, opts.ttel, psfnm, "_sqrtt"*sqrtt, opts.zodi, opts.wfcdark, opts.wfcrnf, opts.wfcrnw, opts.SNRpeaktarg, opts.SNRpeakred, opts.addhost, opts.SNmodel, opts.zlots)

        if not sqrtt:
            fres = open(pltname + ".txt", 'w')

        for i in tqdm.trange(len(phases)):
            for j in tqdm.trange(len(redshifts)):
                plt.subplot(len(phases),len(redshifts),i*len(redshifts) + j + 1)

                SNR5s = []
                SNR10s = []
                SNR20s = []
                SNR50s = []
                SNR90s = []

                for exp_time in exp_times:
                    gal_flambs = make_galaxy_spectrum(redshifts = [redshifts[j]]*nsne)

                    if opts.addhost == 0:
                        gal_flambs = [lambda x:0.] *nsne

                    if prism:
                        SNRs = run_ETC_prism(redshift = redshifts[j], source = source, exp_time = exp_time)
                    else:
                        SNRs = run_ETC(redshift = redshifts[j], phase = phases[i], source = source, exp_time = exp_time, gal_flambs = gal_flambs)

                    SNR5 = scoreatpercentile(SNRs, 5.)
                    SNR50 = scoreatpercentile(SNRs, 50.)
                    SNR20 = scoreatpercentile(SNRs, 20.)
                    SNR10 = scoreatpercentile(SNRs, 10.)
                    SNR90 = scoreatpercentile(SNRs, 90.)

                    print(exp_time, SNR50)

                    SNR5s.append(SNR5)
                    SNR10s.append(SNR10)
                    SNR20s.append(SNR20)
                    SNR50s.append(SNR50)
                    SNR90s.append(SNR90)



                if sqrtt:
                    plt.plot(exp_times, SNR50s/sqrt(exp_times), '.', color = 'k')
                    plt.plot(exp_times, SNR90s/sqrt(exp_times), '.', color = 'b')
                    plt.plot(exp_times, SNR10s/sqrt(exp_times), '.', color = 'r')
                    plt.plot(exp_times, SNR5s/sqrt(exp_times), '.', color = 'm')
                else:
                    plt.plot(exp_times, SNR50s, '.', color = 'k')
                    plt.plot(exp_times, SNR90s, '.', color = 'b')
                    plt.plot(exp_times, SNR10s, '.', color = 'r')
                    plt.plot(exp_times, SNR5s, '.', color = 'm')



                curve50 = fit_curve(exp_times, SNR50s)
                curve90 = fit_curve(exp_times, SNR90s)
                curve5 = fit_curve(exp_times, SNR5s)
                curve10 = fit_curve(exp_times, SNR10s)
                curve20 = fit_curve(exp_times, SNR20s)

                if prism:
                    SNRtargets = [15/2., 20/2., 25/2., 35/2.] + [15./np.sqrt(2.), 20./np.sqrt(2), 25./np.sqrt(2), 35/np.sqrt(2.)]
                    n_dithers_from_SNR = {}
                    for SNRtarget in SNRtargets[:4]:
                        n_dithers_from_SNR[SNRtarget] = 4
                    for SNRtarget in SNRtargets[4:]:
                        n_dithers_from_SNR[SNRtarget] = 2
                else:
                    SNRtargets_one_dither = [opts.SNRpeaktarg*(phases[i] == 0)*np.sqrt((1. + opts.SNRpeakred)/(redshifts[j] + 1.)) +
                                             4.*(phases[i] != 0)*np.sqrt((1. + opts.SNRpeakred)/(redshifts[j] + 1.))]
                    n_dithers_from_SNR = {}
                    for SNRtarget in SNRtargets_one_dither:
                        n_dithers_from_SNR[SNRtarget] = 1

                    SNRtargets = deepcopy(SNRtargets_one_dither)
                    for cadence_to_consider in [1., 2., 2.5, 3., 4., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.]:
                        SNRtargets.append(SNRtargets_one_dither[0]*np.sqrt(cadence_to_consider/5.))
                        n_dithers_from_SNR[SNRtargets[-1]] = 5./cadence_to_consider


                    #SNRtargets = [SNRtargets[0]/np.sqrt(2.)] + SNRtargets
                    #n_dithers_from_SNR[SNRtargets[0]] = 2

                    
                    

                if sqrtt:
                    plt.plot(exp_times, curve90/np.sqrt(exp_times), color = 'b', label = "90th: " + SNR_label(exp_times, curve90, SNRtargets = SNRtargets)[0])
                    plt.plot(exp_times, curve50/np.sqrt(exp_times), color = 'k', label = "50th: " + SNR_label(exp_times, curve50, SNRtargets = SNRtargets)[0])
                    plt.plot(exp_times, curve10/np.sqrt(exp_times), color = 'r', label = "10th: " + SNR_label(exp_times, curve10, SNRtargets = SNRtargets)[0])
                    plt.plot(exp_times, curve5/np.sqrt(exp_times), color = 'm', label = "5th: " + SNR_label(exp_times, curve5, SNRtargets = SNRtargets)[0])
                else:
                    plt.plot(exp_times, curve90, color = 'b', label = "90th: " + SNR_label(exp_times, curve90, SNRtargets = SNRtargets)[0])
                    plt.plot(exp_times, curve50, color = 'k', label = "50th: " + SNR_label(exp_times, curve50, SNRtargets = SNRtargets)[0])
                    plt.plot(exp_times, curve10, color = 'r', label = "10th: " + SNR_label(exp_times, curve10, SNRtargets = SNRtargets)[0])
                    plt.plot(exp_times, curve5, color = 'm', label = "5th: " + SNR_label(exp_times, curve5, SNRtargets = SNRtargets)[0])

                plt.legend(fontsize = 8, loc = 'best')

                if not sqrtt:

                    for thecurve, percentile in ((curve50, 50), (curve20, 20), (curve10, 10), (curve5, 5)):
                        these_results = SNR_label(exp_times, thecurve, SNRtargets = SNRtargets)[1]
                        print("these_results ", these_results)

                        for SNRkey in these_results:
                            n_dithers = n_dithers_from_SNR[SNRkey]
                            print("n_dithers", n_dithers)
                            fres.write("%.1f: filt=%s  TTel=%.1f  redshift=%.1f  phase=%.1f  key=%.2f%s  exp=%.1f\n" % (percentile, opts.filt, opts.ttel, redshifts[j], phases[i], SNRkey,
                                                                                                                        ("x%.3f" % n_dithers)*(n_dithers != 1),
                                                                                                                        these_results[SNRkey]*n_dithers))
                            fres.flush()

                plt.title("S/N @ %i rest-frame, Redshift %.2f" % (phases[i], redshifts[j]))
                plt.xscale('log')
                plt.yscale('log')

        if not sqrtt:
            fres.close()


        plt.savefig(pltname + ".pdf", bbox_inches = 'tight')
        plt.close()
else:
    filts = ["R062", "Z087", "Y106", "J129", "H158", "F184"]
    redshifts = [0.5, 0.9, 1.3, 1.7]

    exp_time = float(sys.argv[4])


    if opts.SNmodel == "salt2_extended":
        source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")
    elif opts.SNmodel == "SALT3.NIR_WAVEEXT":
        source = sncosmo.SALT3Source(modeldir=wfirst_data_path + "/SALT3.NIR_WAVEEXT/")
    else:
        assert 0

    plt.figure(figsize = (6*len(redshifts), 4*6))
    pltname = "StoN_vs_phase_TTel_%.1f_exp=%.1f_%s" % (TTel, exp_time, psfnm)

    for k in range(len(filts)):
        effective_meters2_fl = file_to_fn(wfirst_data_path + "/pixel-level/input/" + filts[k] + ".txt")

        for j in range(len(redshifts)):
            phases = np.arange(-15, 40, 5./(1 + redshifts[j]))
            plt.subplot(len(filts),len(redshifts), k*len(redshifts) + j + 1)

            SNR10s = []
            SNR50s = []
            SNR90s = []

            for i in range(len(phases)):
                print("Redshift", redshifts[j], phases[i])

                gal_flambs = make_galaxy_spectrum(redshifts = [redshifts[j]]*nsne)
                if opts.addhost == 0:
                    gal_flambs = [lambda x:0.] *nsne

                SNRs = run_ETC(redshift = redshifts[j], phase = phases[i], source = source, exp_time = exp_time, gal_flambs = gal_flambs)

                SNR50 = scoreatpercentile(SNRs, 50.)
                SNR20 = scoreatpercentile(SNRs, 20.)
                SNR10 = scoreatpercentile(SNRs, 10.)
                SNR90 = scoreatpercentile(SNRs, 90.)

                SNR10s.append(SNR10)
                SNR50s.append(SNR50)
                SNR90s.append(SNR90)

            SNR10s = np.array(SNR10s)
            SNR50s = np.array(SNR50s)
            SNR90s = np.array(SNR90s)

            stack10 = np.sqrt(np.dot(SNR10s, SNR10s))
            stack50 = np.sqrt(np.dot(SNR50s, SNR50s))
            stack90 = np.sqrt(np.dot(SNR90s, SNR90s))

            plt.plot(phases, SNR10s, '.', color = 'r', label = "10th, stack=%.1f" % stack10)
            plt.plot(phases, SNR50s, '.', color = 'g', label = "50th, stack=%.1f" % stack50)
            plt.plot(phases, SNR90s, '.', color = 'b', label = "90th, stack=%.1f" % stack90)

            plt.title("%s, Redshift %.2f" % (filts[k], redshifts[j]))
            plt.legend(fontsize = 8, loc = 'best')
            plt.xlabel("Phase")
            plt.ylabel("S/N @ Percentile")
            plt.ylim(0, plt.ylim()[1])
    plt.savefig(pltname + ".pdf", bbox_inches = 'tight')
