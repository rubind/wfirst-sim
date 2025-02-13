from numpy import *
import numpy as np
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pixel_level_ETC2 import initialize_PSFs, get_imaging_SN, solve_for_exptime
import sys
from FileRead import readcol, writecol, format_file
import glob
import sncosmo
from DavidsNM import miniLM_new
import tqdm
import os
import sys


def modelfn(P, exp_time, magnitudes):
    """Predict S/N."""
    targ_flux = P[0] * 10.**(-0.4*(magnitudes - 25.)) * exp_time
    back_flux = P[1] * exp_time
    
    nreads = exp_time/3.04

    white_noise_var = P[2]*(nreads - 1)/(nreads + 1)/nreads
    const_noise_var = P[3]

    return targ_flux/sqrt(targ_flux + back_flux + white_noise_var + const_noise_var)

def chi2fn(P, passdata):
    [vals, errs] = passdata[0]

    model = modelfn(P, exp_time, magnitudes)
    return (vals - model)/errs
    
    

def get_mdl(mag, AB_not_ST):
    if AB_not_ST:
        mdl = 10**(-0.4*mag)*0.10884806248/obs_waves**2.
    else:
        mdl = 10**(-0.4*(mag + 21.1)) * ones(len(obs_waves), dtype=float64)
    return mdl


if 0:
    PSF_source = "make_NIC_PSFs"
if 0:
    PSF_source = "Gauss_PSFs"
if 1:
    PSF_source = "WebbPSF_WFC"

AB_not_ST = 1
approximate_PSF = 1
IPC = 0.02
read_noise_white = 16.

wfirst_path = os.environ["WFIRST_SIM_DATA"]


PSFs = initialize_PSFs(pixel_scales = [22], slice_scales = [22], PSF_source = PSF_source)
obs_waves = arange(4000., 25001., 25.)
args = {"waves": obs_waves, "PSFs": PSFs, "gal_flamb": lambda x: 0.0, "TTel": 264., "zodi_fl": "aldering.txt", "source_dir": wfirst_path + "/pixel-level/input/", "read_noise_white": read_noise_white}


magnitudes = arange(21., 29.01, 0.1)
if len(sys.argv) == 1:
    exp_times = 10.**(linspace(log10(20.), log10(2500.), 50))
else:
    exp_times = [float(item) for item in sys.argv[1:]]
    
fexp = open("depth_vs_exptime.txt", 'w')

for effective_meters2_fl, pltcolor in zip(["R062", "Z087", "Y106", "J129", "H158", "F184", "K213"], ['m', 'b', 'c', 'g', 'orange', 'r', 'k']):

    for exp_time in tqdm.tqdm(exp_times):
        vals = []
        errs = []
        
        for magnitude in magnitudes:
            args["mdl"] = get_mdl(magnitude, AB_not_ST)
            tmpvals = []

            for i in range(20):
                o_i = random.random()*22
                o_j = random.random()*22

                ETC_result = get_imaging_SN(redshift = 0.0, exp_time = exp_time, effective_meters2_fl = effective_meters2_fl + ".txt", phase = 0.0, IPC = IPC,
                                            offset_par = o_i, offset_perp = o_j, verbose = False, approximate_PSF = approximate_PSF, **args)
                tmpvals.append(ETC_result["PSF_phot_S/N"])
            vals.append(mean(tmpvals))
            errs.append(  std(tmpvals)/sqrt(float(len(tmpvals)))  )

        vals = array(vals)
        errs = array(errs)

        ifn = interp1d(vals[::-1], magnitudes[::-1], kind = 'linear')
        
        depth = ifn(5.)

        #P, NA, NA = miniLM_new(ministart = [1., 1., 10., 10.], miniscale = [10., 10., 10., 10.], residfn = chi2fn, passdata = [vals, errs])

        #model = modelfn(P, exp_time, magnitudes)
        
        plt.plot(exp_time, depth, '.', color = pltcolor, label = effective_meters2_fl*int(exp_time == exp_times[0]))

        fexp.write("%s %.1f %.2f %.2f %.2f\n" % (effective_meters2_fl, exp_time, depth, depth + 1.25*np.log10(73*0.91), depth + 1.25*np.log10(146*0.91)))
        fexp.flush()
        
            
plt.legend(loc = 'best')
plt.xscale('log')
plt.xlabel("Exposure Time (s)")
plt.ylabel("5-$\sigma$ Depth (AB)")
plt.savefig("5_sigma_vs_time.pdf")
plt.close()

fexp.close()
