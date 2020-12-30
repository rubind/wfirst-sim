from numpy import *
#from FileRead import readcol
import os
import sys
sys.path.append("../pixel-level/")
sys.path.append("../../pixel-level/")
from pixel_level_ETC2 import get_spec_with_err, initialize_PSFs, solve_for_exptime, spectrum_to_matched_resolution, bin_vals_fixed_bins, interpfile
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import commands
from FileRead import readcol, writecol
from string import strip
import glob

wfirst_data_path = os.environ['WFIRST_SIM_DATA']

def get_modjaz_templates():
    f = open("StrippedSNCfAspectra-Modjaz14/CfAspec_StrippedSN_Modjaz_obstable.txt")
    lines = f.read().split('\n')
    f.close()

    selected_sne = []
    selected_dates = []
    for line in lines:
        if len(line) > 42:
            phase = line[35:42]
            try:
                phase = float(phase)
                if phase > -10 and phase < 5:
                    date = line[10:20]
                    selected_dates.append(date.replace(" ", ""))
                    selected_sne.append("sn" + line[3:9].strip())
            except:
                print "Couldn't parse ", phase
                # @ means wrt first spectrum, not LC max.

    print selected_sne
    print selected_dates

    selected_spectra = []
    selected_types = []
    for i in range(len(selected_sne)):
        these_spectra = glob.glob("StrippedSNCfAspectra-Modjaz14/" + selected_sne[i] + "/" + selected_sne[i] + "-" + selected_dates[i] + "*")
        print these_spectra
        selected_spectra += these_spectra
        selected_types += ["CC"]*len(these_spectra)
    return selected_spectra, selected_types


def get_CfA_templates():
    [snparam, NA, tmaxparam] = readcol("cfaspec_snIa/cfasnIa_param.dat", 'aff')
    [snclass, NA, NA, NA, NA, NA, classclass] = readcol("cfaspec_snIa/branchwangclass.dat", 'affffaa')
    [specname, specdate] = readcol("cfaspec_snIa/cfasnIa_mjdspec.dat", 'af')

    spectra2 = []
    sntypes2 = []
    
    print specname[:10]
    print snparam[:10]
    print snclass[:10]
    
    for i in range(len(snparam)):
        found = 0
        for j in range(len(specname)):
            if specname[j].count(snparam[i] + "-"):
                print snparam[i], specname[j]
                found = 1

                if snclass.count(snparam[i]):
                    ind = snclass.index(snparam[i])
                    phase = specdate[j] - tmaxparam[i]
                    
                    if phase > -10 and phase < 5:
                        spectra2.append(glob.glob("cfaspec_snIa/*/" + specname[j])[0])
                        sntypes2.append("Ia-" + classclass[ind])
    return spectra2, sntypes2
                        





spectra, sntypes = get_modjaz_templates()
spectra2, sntypes2 = get_CfA_templates()
spectra += spectra2
sntypes += sntypes2


commands.getoutput("rm -fr simulated_observations")
for sntype in set(sntypes):
    commands.getoutput("mkdir -p simulated_observations/" + sntype)


#[spectra, types] = readcol("selected_spectra.txt", 'aa')
PSFs = initialize_PSFs(pixel_scales = [15, 30], slice_scales = [30, 30], PSF_source = "WebbPSF", path = wfirst_data_path + "/pixel-level/")

args = {"pixel_scale": 0.05, "slice_scale": 0.15, "source_dir": wfirst_data_path + "/pixel-level/input/", "IFURfl": "IFU_R_160720.txt", "min_wave": 4200.,
        "max_wave": 21000., "PSFs": PSFs, "phase": 0, "redshift": 1.0}


first_ETC_run = get_spec_with_err(exp_time = 100., **args)
args["waves"] = first_ETC_run["obs_waves"]

rest_waves = args["waves"]/(1 + args["redshift"])

for spectrum, sntype in zip(spectra, sntypes):
    print spectrum

    [orig_waves, orig_fluxes] = readcol(spectrum, 'ff')
    orig_fluxes /= max(orig_fluxes)

    # First, we just get the correct normalization

    spec = interp1d(orig_waves, orig_fluxes, kind = 'linear', fill_value = 0, bounds_error = False)(rest_waves)

    inds = where((rest_waves > 5000.)*(rest_waves < 6000.))
    orig_fluxes *= dot(first_ETC_run["f_lamb_SN"][inds], first_ETC_run["f_lamb_SN"][inds])/dot(spec[inds], first_ETC_run["f_lamb_SN"][inds])
    
    matched_resolution = spectrum_to_matched_resolution(orig_waves*(1 + args["redshift"]), orig_fluxes, min_wave = args["min_wave"], max_wave = args["max_wave"],
                                                        IFUR = interpfile(wfirst_data_path + "/pixel-level/input/IFU_R_160720.txt"), pixel_scale = 0.075)


    args["mdl"] = matched_resolution
    exp_time = solve_for_exptime(S_to_N = 23.2, key1= "rest_frame_band_S/N", key2= (5000, 6000), **args)

    # Now, we can run the ETC


    ETC_run = get_spec_with_err(exp_time = exp_time, **args)
    


    flux = ETC_run["f_lamb_SN"]
    fluxerr = ETC_run["f_lamb_SN"]/ETC_run["spec_S/N"]

    good_inds = where((1 - isnan(flux))*(1 - isnan(fluxerr))*(1 - isinf(fluxerr))*(rest_waves > 3750.)*(rest_waves < 7600.))


    writecol("simulated_observations/" + sntype + "/" + spectrum.split("/")[-1].split(".")[0] + ".txt", [args["waves"][good_inds], flux[good_inds] + random.normal(size = len(fluxerr[good_inds]))*fluxerr[good_inds], fluxerr[good_inds]])
