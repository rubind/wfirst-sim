from copy import deepcopy
from numpy import *
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pixel_level_ETC2 import initialize_PSFs, get_imaging_SN, solve_for_exptime, get_spec_with_err
import sys
from FileRead import readcol, writecol, format_file
import commands
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

def file_to_fn(fl, col = 1):
    vals = loadtxt(fl)
    x = vals[:,0]
    y = vals[:,col]

    return interp1d(x, y, kind = 'linear')


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
    [MB, x1, color, mass] = [item[0] for item in MB, x1, color, mass]

    MB += random.normal()*0.1 + 0.055*redshift*random.normal() + 5/log(10.)*(0.001/redshift)*random.normal()

    sncosmo_model = get_SNCosmo_model(redshift = redshift,
                                      x1 = x1,
                                      c = color,
                                      MV = MB - color,
                                      daymax = daymax, source = source)


    return MB - color, x1, color, mass, sncosmo_model

def run_ETC(redshift, phase, source, exp_time):
    SNRs = []

    for i in range(nsne):
        NA, NA, NA, NA, sncosmo_model = realize_SN(redshift = redshift, daymax = 0, source = source)

        f_lamb_SN = sncosmo_model.flux(phase*(1. + redshift), WFI_args["waves"])
        
        ETC_result = get_imaging_SN(redshift = 0, exp_time = exp_time, gal_flamb = gal_flambs[i],
                                    effective_meters2_fl = effective_meters2_fl, phase = 0, mdl = f_lamb_SN,
                                    offset_par = random.random()*22, offset_perp = random.random()*22, **WFI_args)
        SNRs.append(ETC_result["PSF_phot_S/N"])

    return SNRs

def make_SN_curve(P, exp_time):
    return exp_time/sqrt(P[0] + P[1]/exp_time + P[2]*exp_time)

def chi2fn(P, thedata):
    exp_times, SNs = thedata[0]
    
    if any(P < 0):
        return 1e100

    model = make_SN_curve(P, exp_times)
    resid = (SNs - model)/SNs
    return dot(resid, resid)

def fit_curve(exp_times, SNs):
    bestF = 1e100

    for ministart in ([1., 1., 1.],
                      [10., 10., 10.],
                      [0.2, 0.2, 0.2],
                      [50., 50., 50.]):
        

        P, F, NA = miniNM_new(ministart = ministart, miniscale = [0.1, 0.1, 0.1], chi2fn = chi2fn, passdata = [exp_times, SNs])
    
        if F < bestF:
            bestF = F
            bestP = P

    return make_SN_curve(bestP, exp_times)

def SNR_label(exp_times, curve, SNRtargets):
    label = ""
    vals = {}

    for SNRtarget in SNRtargets:
        if curve.min() < SNRtarget and curve.max() > SNRtarget:
            ifn = interp1d(curve, exp_times, kind = 'linear')
            label += "S/N %.1f: %.2g s " % (SNRtarget, ifn(SNRtarget))
            vals[SNRtarget] = ifn(SNRtarget)
    return label, vals

print "Reading PSFs..."
TTel = float(sys.argv[2])
psfnm = "WebbPSF_WFC"

PSFs = initialize_PSFs(pixel_scales = [22], slice_scales = [22], PSF_source = psfnm)
WFI_args = {"PSFs": PSFs, "source_dir": wfirst_data_path + "/pixel-level/input/",
            "pixel_scale": 0.11, "dark_current": 0.015,
            "IPC": 0.02,
            "TTel": TTel,
            "zodi_fl": file_to_fn(wfirst_data_path + "/pixel-level/input/aldering.txt"),
            "bad_pixel_rate": 0.01,
            "waves": arange(6000., 23000.1, 25.)}


redshifts = [0.4, 0.8, 1.7, 2.0]
phases = [-10, 0]
nsne = 40
filt = sys.argv[1]

exp_times = 10**arange(1.0, 4.01, 0.05)
source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")
effective_meters2_fl = file_to_fn(wfirst_data_path + "/pixel-level/input/" + filt + ".txt")

for sqrtt in [0, 1]:
    plt.figure(figsize = (6*len(redshifts), 4*len(phases)))
    pltname = "StoN_vs_exptime_%s_TTel_%.1f_%s%s" % (filt, TTel, psfnm, "_sqrtt"*sqrtt)

    if not sqrtt:
        fres = open(pltname + ".txt", 'w')

    for i in range(len(phases)):
        for j in range(len(redshifts)):
            plt.subplot(len(phases),len(redshifts),i*len(redshifts) + j + 1)

            SNR10s = []
            SNR50s = []
            SNR90s = []
            
            for exp_time in exp_times:
                gal_flambs = make_galaxy_spectrum(redshifts = [redshifts[j]]*nsne)

                SNRs = run_ETC(redshift = redshifts[j], phase = phases[i], source = source, exp_time = exp_time)

                SNR50 = scoreatpercentile(SNRs, 50.)
                SNR10 = scoreatpercentile(SNRs, 10.)
                SNR90 = scoreatpercentile(SNRs, 90.)

                print exp_time, SNR50

                SNR10s.append(SNR10)
                SNR50s.append(SNR50)
                SNR90s.append(SNR90)



            if sqrtt:
                plt.plot(exp_times, SNR50s/sqrt(exp_times), '.', color = 'k')
                plt.plot(exp_times, SNR90s/sqrt(exp_times), '.', color = 'b')
                plt.plot(exp_times, SNR10s/sqrt(exp_times), '.', color = 'r')
            else:
                plt.plot(exp_times, SNR50s, '.', color = 'k')
                plt.plot(exp_times, SNR90s, '.', color = 'b')
                plt.plot(exp_times, SNR10s, '.', color = 'r')



            curve50 = fit_curve(exp_times, SNR50s)
            curve90 = fit_curve(exp_times, SNR90s)
            curve10 = fit_curve(exp_times, SNR10s)
                
            SNRtargets = [10.*(phases[i] == 0) + 4.*(phases[i] != 0)]
            SNRtargets = [SNRtargets[0]/sqrt(2.)] + SNRtargets # Only considering two dithers below, so don't change this

            if sqrtt:
                plt.plot(exp_times, curve90/sqrt(exp_times), color = 'b', label = "90th: " + SNR_label(exp_times, curve90, SNRtargets = SNRtargets)[0])
                plt.plot(exp_times, curve50/sqrt(exp_times), color = 'k', label = "50th: " + SNR_label(exp_times, curve50, SNRtargets = SNRtargets)[0])
                plt.plot(exp_times, curve10/sqrt(exp_times), color = 'r', label = "10th: " + SNR_label(exp_times, curve10, SNRtargets = SNRtargets)[0])
            else:
                plt.plot(exp_times, curve90, color = 'b', label = "90th: " + SNR_label(exp_times, curve90, SNRtargets = SNRtargets)[0])
                plt.plot(exp_times, curve50, color = 'k', label = "50th: " + SNR_label(exp_times, curve50, SNRtargets = SNRtargets)[0])
                plt.plot(exp_times, curve10, color = 'r', label = "10th: " + SNR_label(exp_times, curve10, SNRtargets = SNRtargets)[0])
            
            plt.legend(fontsize = 8, loc = 'best')
            
            if not sqrtt:

                for thecurve, percentile in ((curve50, 50), (curve10, 10)):
                    these_results = SNR_label(exp_times, thecurve, SNRtargets = SNRtargets)[1]

                    for SNRkey in these_results:
                        isminexptime = SNRkey == min(these_results.keys())
                        
                        fres.write("%.1f: filt=%s  TTel=%.1f  redshift=%.1f  phase=%.1f  key=%.2f%s  exp=%.1f\n" % (percentile, filt, TTel, redshifts[j], phases[i], SNRkey,
                                                                                                                    "x2"*isminexptime, these_results[SNRkey]*(1 + isminexptime)))
                        fres.flush()

            plt.title("S/N @ %i rest-frame, Redshift %.2f" % (phases[i], redshifts[j]))
            plt.xscale('log')
            plt.yscale('log')

    if not sqrtt:
        fres.close()
        

    plt.savefig(pltname + ".pdf", bbox_inches = 'tight')
    plt.close()

