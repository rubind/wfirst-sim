from copy import deepcopy
from numpy import *
from DavidsNM import miniNM_new, save_img
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.environ['WFIRST'] + "/scripts/pixel-level/")
from pixel_level_ETC2 import get_spec_with_err, initialize_PSFs, solve_for_exptime
from scipy.interpolate import interp1d
from SimpleFisher import run_FoM, get_Jacobian_NSNe1, get_params
import argparse
from astropy.cosmology import FlatLambdaCDM

def volume_of_z(z):
    return interp1d(
        [0.0, 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8],
        [0.0, 0.33, 2.48, 7.74, 16.9, 30.5, 48.5, 70.9, 97.5, 128., 162., 199., 239., 281., 325., 371., 419., 468., 511.],
        kind = 'linear')(z)

def sn_rates(z):
    # 1e-4 /year/Mpc^3/h70^3.
    margin = 0.9

    return interp1d([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5],
                    [0.218295627331*margin,  0.249056423859*margin,  0.287104418451*margin,  0.329567530646*margin,  0.376198210965*margin,  0.426692279263*margin,  0.479950657435*margin,  0.534831964795*margin,  0.589178179388*margin,  0.640535383579*margin,  0.686037612282*margin,  0.723178391194*margin,  0.749490710152*margin,  0.763949149987*margin,  0.765811346301*margin,  0.756043123009*margin,  0.736584324802*margin,  0.709373826398*margin,  0.676748093973*margin,  0.640659438231*margin,  0.602657314567*margin,  0.564327310361*margin,  0.526540587149*margin,  0.489816696065*margin,  0.454676544907*margin,  0.421393249604*margin],
                    kind = 'linear')(z)

def sne_in_sqdeg_per_twoobsyear(z, dz = 0.1):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    sne_in_full_sky_per_rest_year = cosmo.comoving_volume(z + dz/2.).value - cosmo.comoving_volume(z - dz/2.).value
    sne_in_full_sky_per_rest_year *= sn_rates(z)*1.e-4
    sne_in_full_sky_per_twoobs_year = 2.*sne_in_full_sky_per_rest_year/(1. + z)
    sne_in_sqdeg_per_twoobs_year = sne_in_full_sky_per_twoobs_year/41252.96
    return sne_in_sqdeg_per_twoobs_year

def search_time(z):

    return interp1d([0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0],
                    [21.6 + 2.*opts.slew, 21.6 + 2.*opts.slew, 37.5 + 2.*opts.slew, 52.9 + 2.*opts.slew, 65.5 + 2.*opts.slew,
                     81.2 + 2.*opts.slew, 94.1 + 2.*opts.slew, 112.5 + 2.*opts.slew, 129.6 + 2.*opts.slew, 148.9 + 2.*opts.slew,
                     172.5 + 2.*opts.slew, 205.6 + 2.*opts.slew, 236.1 + 2.*opts.slew, 277.4 + 2.*opts.slew, 307.2 + 2.*opts.slew,
                     347.6 + 2.*opts.slew, 401.8 + 2.*opts.slew, 454.8 + 2.*opts.slew, 535.0 + 2.*opts.slew, 629.3 + 2.*opts.slew], kind = 'linear')(z)


def supernova_survey_time(redshifts, sn_counts, verbose = False):
    sne_found = zeros(len(redshifts), dtype=float64)
    exp_time = 0

    if verbose:
        plt.figure()
        plt.plot(redshifts, sn_counts, 'o')

    for i in range(len(redshifts))[::-1]:
        sne_to_find = sn_counts[i] - sne_found[i]

        if sne_to_find > 0:
            sqdeg = sne_to_find/sne_in_sqdeg_per_twoobsyear(redshifts[i])
            pointings = sqdeg / 0.28
            exp_time += pointings*search_time(redshifts[i])*146. # About 140 search visits
            sne_found[i] += sne_in_sqdeg_per_twoobsyear(redshifts[i])*sqdeg

            if verbose:
                print "Surveying at z=" + str(redshifts[i]) + " sqdeg=" + str(sqdeg)

            for j in range(0,i):
                sne_found[j] += sne_in_sqdeg_per_twoobsyear(redshifts[j])*sqdeg
            
            if verbose:
                plt.plot(redshifts, sne_found, label = "%.2f  so far: %.2g s %.2g d" % (redshifts[i], exp_time, exp_time/86400.))

    if verbose:
        plt.legend(loc = 'best')
        plt.savefig("sn_survey.pdf")
        plt.close()
    #print "HACK!!!!!"*100
    return exp_time#*2.5

"""
class tmpclass:
    def __init__(self):
        self.slew = 70.

opts = tmpclass()
supernova_survey_time(arange(0.15, 1.66, 0.1), [float(item) for item in "60.77459864  159.32710831  291.06756882  451.99380599  618.66594294  260.02681612  177.87616905  177.62276607  187.9598534   206.62014954  208.51613857  201.82254882  193.18932657  167.3886341   147.07147408  124.02081812".split(None)], verbose = True)
stop_here
"""

def make_paramfile(scaled_guess):
    f = open("paramfile_wrap.txt")
    orig_lines = f.read()
    f.close()

    lines = orig_lines.replace("NNNNN", str([800.] + list(scaled_guess)))
    lines = lines.replace("FFFFF", str([opts.FoM] + opts.bins))
    lines = lines.replace("EEEEE", str(at_max_exp_times))
    f = open("paramfile_tmp.txt", 'w')
    f.write(lines)
    f.close()


def chi2fn(new_guess, NA):
    eval_time = dot(exp_times, new_guess)
    survey_time = supernova_survey_time(redshifts, new_guess)*(1 - find_from_ground)
    eval_time += survey_time

    scaled_guess = new_guess*total_time/eval_time

    print "scaled ", list(scaled_guess)
    
    these_params = deepcopy(params)
    these_params["NSNe"] = concatenate(([800.], scaled_guess))

    FoM = run_FoM(jacobian_NSNe1 = jacobian_NSNe1, PSFs = PSFs, params = these_params)
    print "FoM resulting from scaled: ", FoM

    return -FoM


parser = argparse.ArgumentParser()
parser.add_argument('-slew', help='Slew time (s)', type = float)
parser.add_argument('-maxz', help='Maximum Redshift', type = float, default = 1.8)
parser.add_argument('-FoM', help='FoM type', type = str)
parser.add_argument('-bins', help='Redshift bins, starts with 0', type = float, default = [], nargs = '*')
opts = parser.parse_args()


find_from_ground = 0


PSFs = initialize_PSFs(pixel_scales = [10], slice_scales = [30], PSF_source = "WebbPSF")
redshifts = arange(0.15, opts.maxz - 0.04, 0.1)
exp_times = []
at_max_exp_times = []

for redshift in redshifts:
    exp_time = solve_for_exptime(10.*sqrt(15.), redshift, PSFs, key1 = "rest_frame_band_S/N", key2 = (5000, 6000),
                                 pixel_scale = 0.05, slice_scale = 0.15,
                                 source_dir = os.environ["WFIRST_SIM_DATA"] + "/pixel-level/input/",
                                 IFURfl = "IFU_R_160720.txt", min_wave = 4200.)*3
    at_max_exp_times.append(exp_time/3.)

    exp_time += opts.slew*6. # 6 visits; 1+1 + 4-point ref
    exp_times.append(exp_time)
    if 0:
        exp_times.append(0.)
        print "HACK!!!!!"*100

exp_times = array(exp_times)
print "exp_times ", exp_times

plt.plot(redshifts, exp_times)

plt.savefig("exptime_vs_z.pdf")
plt.close()


make_paramfile(scaled_guess = 100*ones(len(redshifts), dtype=float64))
params = get_params("paramfile_tmp.txt", PSFs)
jacobian_NSNe1 = get_Jacobian_NSNe1(params)

initial_guess = exp(-redshifts/2.)*10. * (1 + random.random(size = len(redshifts)))

total_time = 15778463.04


P, F, NA = miniNM_new(ministart = initial_guess, miniscale = initial_guess/3., chi2fn = chi2fn, passdata = None, verbose = True, inlimit = lambda x: all(x >= 0), maxruncount = 100, maxiter = 500)


