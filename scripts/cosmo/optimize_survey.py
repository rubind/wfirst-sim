from numpy import *
from DavidsNM import miniNM_new
import matplotlib.pyplot as plt
import sys
sys.path.append("../pixel-level/")
import os
from pixel_level_ETC2 import get_spec_with_err, initialize_PSFs, solve_for_exptime
from scipy.interpolate import interp1d
from SimpleFisher import run_FoM

def volume_of_z(z):
    return interp1d(
        [0.0, 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7],
        [0.0, 0.33, 2.48, 7.74, 16.9, 30.5, 48.5, 70.9, 97.5, 128., 162., 199., 239., 281., 325., 371., 419., 468.],
        kind = 'linear')(z)

def sn_rates(z):
    # 1e-4 /year/Mpc^3/h70^3.
    return interp1d([0., 1., 2.], [0.25, 0.75, 0.5], kind = 'linear')(z)

def sne_in_sqdeg_per_twoobsyear(z, dz = 0.1):
    sne_in_full_sky_per_rest_year = volume_of_z(z + dz/2.) - volume_of_z(z - dz/2.)
    sne_in_full_sky_per_rest_year *= 1.e9*sn_rates(z)*1.e-4
    sne_in_full_sky_per_twoobs_year = 2.*sne_in_full_sky_per_rest_year/(1. + z)
    sne_in_sqdeg_per_twoobs_year = sne_in_full_sky_per_twoobs_year/41252.96
    return sne_in_sqdeg_per_twoobs_year

def search_time(z):

    return interp1d([0.,
                     0.8,
                     1.7],
                    [18.*2 + 70*2.,
                     67.*2 + 70*2.,
                     265.*2 + 70.*2], kind = 'linear')(z)


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
            exp_time += pointings*search_time(redshifts[i])*146. # About 140 search orbits
            sne_found[i] += sne_in_sqdeg_per_twoobsyear(redshifts[i])*sqdeg

            if verbose:
                print "Surveying at z=" + str(redshifts[i]) + " sqdeg=" + str(sqdeg)

            for j in range(0,i):
                sne_found[j] += sne_in_sqdeg_per_twoobsyear(redshifts[j])*sqdeg
            
            if verbose:
                plt.plot(redshifts, sne_found, label = "%.2f  %.2g s %.2g d" % (redshifts[i], exp_time, exp_time/86400.))

    if verbose:
        plt.legend(loc = 'best')
        plt.savefig("sn_survey.pdf")
        plt.close()
    return exp_time/5.

#supernova_survey_time(arange(0.15, 1.66, 0.1), [float(item) for item in "60.77459864  159.32710831  291.06756882  451.99380599  618.66594294  260.02681612  177.87616905  177.62276607  187.9598534   206.62014954  208.51613857  201.82254882  193.18932657  167.3886341   147.07147408  124.02081812".split(None)], verbose = True)
#supernova_survey_time(arange(0.15, 1.66, 0.1), [100]*16, verbose = True)

def chi2fn(new_guess, NA):
    eval_time = dot(exp_times, new_guess)
    survey_time = supernova_survey_time(redshifts, new_guess)*(1 - find_from_ground)
    eval_time += survey_time

    scaled_guess = new_guess*total_time/eval_time

    print "scaled ", list(scaled_guess)
    
    lines = orig_lines.replace("NNNNN", str([800.] + list(scaled_guess)))
    lines = orig_lines.replace("FFFFF", '["' + str(sys.argv[3]) + '"]')
    f = open("paramfile_tmp.txt", 'w')
    f.write(lines)
    f.close()

    FoM = run_FoM("paramfile_tmp.txt", PSFs = PSFs)
    print "FoM resulting from scaled: ", FoM

    return -FoM



find_from_ground = 0

f = open("paramfile_wrap.txt")
orig_lines = f.read()
f.close()

PSFs = initialize_PSFs(pixel_scales = [10], slice_scales = [30], PSF_source = "WebbPSF", path = "../pixel-level/")
redshifts = arange(0.15, 1.66, 0.1)
exp_times = []
for redshift in redshifts:
    exp_time = solve_for_exptime(10., redshift, PSFs, key1 = "obs_frame", key2 = (10200, 12850), pixel_scale = 0.05, slice_scale = 0.15,
                                 source_dir = "../pixel-level/input/", IFURfl = "IFU_R_160720.txt", min_wave = 4200.)*3
    exp_time += 70.*6 # 6 visits; 1+1 + 4-point ref
    exp_times.append(exp_time)

exp_times = array(exp_times)
print "exp_times ", exp_times

plt.plot(redshifts, exp_times)

plt.savefig("exptime_vs_z.pdf")
plt.close()

initial_guess = ones(16.) #array([69, 208, 402, 223, 327, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136.])

total_time = 15778463.04


P, F, NA = miniNM_new(ministart = initial_guess, miniscale = initial_guess/3., chi2fn = chi2fn, passdata = None, verbose = True, inlimit = lambda x: all(x > 0))


