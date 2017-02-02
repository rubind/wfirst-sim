from numpy import *
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import commands
import sys


def get_filters(include_ground, is_deep, spec_SNR, TTel):
    if TTel < 270:
        reddest_filt = "F184"
        reddest_exp = 973.4
    else:
        reddest_filt = "F184"
        reddest_exp = 2312.9

    if include_ground:
        if is_deep:
            # Deep tier, with ground-based data: YJ + HF parallels (or YJHF if no spectroscopy)
            n_WFC_filters = 2 + 2*(spec_SNR == None)
            filters = ["g", "r", "i", "z", "Y"] + ["Y106", "J129"] + ["H158", reddest_filt]*(spec_SNR == None)
            exp_times = [600, 600, 600, 1800, 1800] + [243.3, 268.8] + [361.6, reddest_exp]*(spec_SNR == None)
        else:
            # Shallow tier, with ground-based data: nothing (H if no spectroscopy)
            n_WFC_filters = (spec_SNR == None)
            filters = ["g", "r", "i", "z", "Y"] + ["H158"]*(spec_SNR == None)
            exp_times = [600, 600, 600, 1800, 1800] + [137.6]*(spec_SNR == None)
    else:
        if is_deep:
            # Deep tier, without ground-based data: YJ + HF parallels (or YJHF if no spectroscopy)
            n_WFC_filters = 2 + 2*(spec_SNR == None)
            filters = ["Y106", "J129"] + ["H158", reddest_filt]*(spec_SNR == None)
            exp_times = [243.3, 268.8] + [361.6, reddest_exp]*(spec_SNR == None)
        else:
            # Shallow tier, without ground-based data: ZYJ + H parallels (or ZYJH if no spectroscopy)
            n_WFC_filters = 3 + 1*(spec_SNR == None)
            filters = ["Z087", "Y106", "J129"] + ["H158"]*(spec_SNR == None)
            exp_times = [50.7, 55.0, 86.8] + [137.6]*(spec_SNR == None)

    assert len(filters) == len(exp_times)
    return int(n_WFC_filters), filters, exp_times, exp_times[-int(n_WFC_filters):], reddest_filt



def make_param_file(flname, include_ground, spec_SNR, deep_square_degrees, shallow_square_degrees, frac_deep, max_deep_z, slope_to_max_z, deep_scale, spectra_depth_and_phase, big_survey, h1rg, survey_time, TTel):
    working_dir = flname[:flname.rfind("/")]
    pwd = commands.getoutput("pwd")
    just_flname = flname.split("/")[-1]
    assert just_flname != flname

    commands.getoutput("mkdir -p " + working_dir)

    f = open(working_dir + "/survey.sh", 'w')
    f.write("#!/bin/bash -l\n")
    f.write("#SBATCH --partition=shared\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --time=24:00:00\n")
    f.write("#SBATCH --job-name=survey\n")
    f.write("#SBATCH --mem=%i00\n" % (25 + 40*big_survey + 40*(shallow_square_degrees > 40) ))
    f.write("#SBATCH -L SCRATCH\n")
    hostname = commands.getoutput("hostname")
    if hostname.count("cori"):
        f.write("#SBATCH -C haswell\n")
    f.write("module load python/2.7-anaconda\n")
    f.write("cd " + pwd + "/" + working_dir + "\n")
    wfirst = commands.getoutput("echo $WFIRST")
    f.write("srun -n 1 -c 1 python " + wfirst + "/scripts/stan_cosmo/STEP1_simulate_survey.py " + just_flname + " " + just_flname.replace(".csv", ".txt").replace("paramfile", "pickle") + " > log1.txt\n")
    f.close()

    f = open(flname, 'w')
    if big_survey == 0:
        f.write("total_survey_time,%.2f\n" % survey_time)
        f.write("survey_duration,2.5\n") # Separate from how long the survey actually takes. This is the maximum!
    elif big_survey == 1:
        f.write("total_survey_time,%.2f\n" % (survey_time*2.))
        f.write("survey_duration,5\n")
    else:
        f.write("total_survey_time,%.3f\n" % (survey_time*3.))
        f.write("survey_duration,7.5\n")

    f.write("maximum_trigger_fraction,0.9\n")					
    f.write("adjust_each_SN_exp_time,FALSE\n")					
    f.write("normalization_wavelength_range,5000,6000\n")
    
    if spec_SNR != None:
        f.write("shallow_SNR," + spec_SNR[0] + "\n")
        f.write("medium_SNR," + spec_SNR[1] + "\n")
        f.write("deep_SNR," + spec_SNR[2] + "\n")
        f.write("reference_SNR," + spec_SNR[3] + "\n")
        f.write("spectra_depth_and_phase," + ",".join(spectra_depth_and_phase) + "\n")
    else:
        f.write("shallow_SNR,10\n")
        f.write("medium_SNR,10\n")
        f.write("deep_SNR,10\n")
        f.write("reference_SNR,10\n")
        f.write("spectra_depth_and_phase,,\n")
        
    f.write("targeted_parallels,TRUE\n")
    f.write("number_of_reference_dithers,4\n")
    f.write("SN_rates,SN_rates.txt\n")					
    f.write("grizY_30s_ground_depths,24.47,24.16,23.4,22.23,21.57\n")
    f.write("__________________________________________\n")
    f.write("slew_table,slewtable_I32000_t0.40_vmax0.12_Tr2.0.txt\n")					
    f.write("zodiacal_background,aldering.txt\n")		
    f.write("telescope_temperature,%.0f\n" % TTel)	
    f.write("PSFs,WebbPSF\n")
    f.write("IFU_read_noise_floor,4.0\n")
    f.write("interpixel_capacitance,0.02\n")	
    f.write("WFI_dark_current,0.015\n")	
    f.write("IFU_min_wave,4200\n")
    f.write("IFU_max_wave,21000\n")
    f.write("IFU_effective_area,IFU_effective_area_160720.txt\n")
    #f.write("IFU_resolution,IFU_R_Content.txt\n")	
    f.write("IFU_resolution,IFU_R_160720.txt\n")
    f.write("IFU_pixel_scale,%s\n" % (h1rg*"0.075" + (1 - h1rg)*"0.05"))
    f.write("IFU_slice_in_pixels,%s\n" % (h1rg*"2" + (1 - h1rg)*"3"))
    f.write("IFU_dark_current,%s\n" % (h1rg*"0.01" + (1 - h1rg)*"0.003"))
    f.write("bad_pixel_rate,0.01\n")

    if include_ground:
        survey_tiers = ["Wide1", "Wide2", "Wide3"] + ["Deep"]*(deep_square_degrees > 0)
    else:
        survey_tiers = ["Wide"] + ["Deep"]*(deep_square_degrees > 0)
    

    for tier in survey_tiers:
        is_deep = tier == "Deep"
        n_shallow = 1. + 2*include_ground

        f.write("__________________________________________\n")
        
        f.write("tier_name," + tier + "\n")
        f.write("tier_fraction_time,%f\n" % (is_deep*frac_deep + (1 - is_deep)*(1 - frac_deep)/n_shallow))
        f.write("square_degrees,%.2f\n" % (deep_square_degrees*is_deep + (1 - is_deep)*shallow_square_degrees/n_shallow))


        n_WFC_filters, filters, exp_times_per_dither, NA, reddest_filt = get_filters(include_ground, is_deep, spec_SNR, TTel)

        f.write("parallel_filters,%s\n" % ("H158" + ("," + reddest_filt)*is_deep))




        dithers_per_filter = [1]*(5*include_ground) + [1]*n_WFC_filters
        dithers_per_filter = [str(item) for item in dithers_per_filter]
        exp_times_per_dither = [str(item) for item in exp_times_per_dither]


        f.write("filters," + ",".join(filters) + "\n")
        f.write("exp_times_per_dither," + ",".join(exp_times_per_dither) + "\n")
        f.write("cadence,5\n")
        f.write("dithers_per_filter," + ",".join(dithers_per_filter) + "\n")

        if not is_deep:
            if include_ground:
                trigger_redshifts = [0, 0.5, 0.65, 0.8, 0.81, 3]
                trigger_fraction = [1, 1, 0.6, 0.3, 0, 0]
            else:
                trigger_redshifts = [0, 0.5, 0.65, 0.8, 0.81, 3]
                trigger_fraction = [1, 1, 1, 1, 0, 0]

        else:
            estimated_fraction_to_target = clip(deep_scale/deep_square_degrees, 0, 1)
            if max_deep_z < 1.6:
                estimated_fraction_to_target *= 1.5

            trigger_redshifts = [0, max_deep_z*0.8,  max_deep_z, max_deep_z + 0.01, 3]
            trigger_fraction = [estimated_fraction_to_target, estimated_fraction_to_target, estimated_fraction_to_target*(1. - slope_to_max_z*0.4), 0, 0]
        trigger_redshifts = ["%.2f" % item for item in trigger_redshifts]
        trigger_fraction = ["%.2f" % item for item in trigger_fraction]

        f.write("trigger_redshifts," + ",".join(trigger_redshifts) + "\n")
        f.write("trigger_fraction," + ",".join(trigger_fraction) + "\n")

        
    f.close()
    print commands.getoutput("sbatch " + working_dir + "/survey.sh")


def lnspace(start, stop, samples):
    return list(around(exp(linspace(log(start), log(stop), samples))))

def sqspace(start, stop, samples, decimals = 0):
    return list(around(linspace(sqrt(start), sqrt(stop), samples)**2., decimals = decimals))

def get_combination():
    params = {}

    params["include_ground"] = 0 #random.choice([0, 1], size = 1)[0]
    #params["spec_SNR"] = random.choice([("12", "20.6", "34.3", "40"), ("13.5", "23.2", "38.6", "45"), None], size = 1, p = [0.425, 0.425, 0.15])[0]
    params["spec_SNR"] = [3.5*sqrt(15.), 6*sqrt(15.), 10.*sqrt(15), sqrt(136.*15.)]
    params["spec_SNR"] = ["%.1f" % item for item in params["spec_SNR"]]

    params["deep_square_degrees"] = 5.#random.choice(  concatenate(([0], sqspace(3, 10, 3))) , size = 1)[0]
    if params["deep_square_degrees"] > 0:
        params["deep_scale"] = random.choice(  [1.0] , size = 1)[0]
    else:
        params["deep_scale"] = 0

    if params["deep_square_degrees"] > 0:
        if params["include_ground"]:
            params["shallow_square_degrees"] = 9.6*3
        else:
            params["shallow_square_degrees"] = random.choice(  sqspace(10, 30, 5) , size = 1)[0]
    else:
        # Need about 100 square degrees
        params["shallow_square_degrees"] = 96.

    if params["deep_square_degrees"] > 0:
        params["frac_deep"] = random.choice([0.6, 0.7, 0.8], size = 1)[0]
    else:
        params["frac_deep"] = 0

    params["max_deep_z"] = random.choice([1.1, 1.4, 1.7, 2.0], size = 1)[0]
    params["slope_to_max_z"] = random.choice([0, 1], size = 1)[0]

    params["spectra_depth_and_phase"] = ("shallow  -10", "medium -1", "deep 1", "shallow 5", "shallow 10")
    
    
    
    params["h1rg"] = 1#random.choice([0, 1], size = 1)[0]
    params["TTel"] = random.choice([260], size = 1)[0]

    params["survey_time"] = 0.6#random.choice([0.5, 0.6], size = 1)[0]




    # Let's do some basic checking: do we go over 12 hours of cadenced imaging with ~40 seconds dither time and two filters * two dithers (four filters without spectra)?

    NA, NA, NA, exp_times_shallow, NA = get_filters(include_ground = params["include_ground"], is_deep = 0, spec_SNR = params["spec_SNR"], TTel = params["TTel"])

    NA, NA, NA, exp_times_deep, NA = get_filters(include_ground = params["include_ground"], is_deep = 1, spec_SNR = params["spec_SNR"], TTel = params["TTel"])
    
    number_of_deep_pointings = around(params["deep_square_degrees"]/0.28)

    if params["include_ground"]:
        number_of_shallow_pointings = 0
    else:
        number_of_shallow_pointings = around(params["shallow_square_degrees"]/0.28)

    total_exp_time = number_of_deep_pointings*sum(array(exp_times_deep) + 62.)  + number_of_shallow_pointings*sum(array(exp_times_shallow) + 62.)

    if params["spec_SNR"] != None:
        if total_exp_time > 14*3600:
            print "Regenerating..."
            return get_combination()
    else:
        if total_exp_time < 27*3600 or total_exp_time > 33*3600:
            print "Regenerating..."
            return get_combination()


    return params


nsurveys = 200

combs = [get_combination() for i in range(nsurveys)]

print len(combs)


params_as_array = []


for item in combs:
    snr = item["spec_SNR"]
    if snr == None:
        snr = 0.
    else:
        snr = float(snr[0])

    fd = item["frac_deep"]
    if fd == None:
        fd = 0

    params_as_array.append([item["max_deep_z"],
                            item["shallow_square_degrees"],
                            item["h1rg"],
                            snr,
                            item["deep_square_degrees"],
                            item["include_ground"],
                            item["slope_to_max_z"]])

    


params_as_array = array(params_as_array)

print params_as_array.shape

cmat = corrcoef(transpose(params_as_array))

print cmat.shape

print combs[0]

from DavidsNM import save_img
save_img(cmat, "cmat.fits")




commands.getoutput("rm -f slurm-*.out")
commands.getoutput("rm -fr " + sys.argv[1])


surveys_written = 0



print "Ready to start..."

for i in range(nsurveys):
    make_param_file(sys.argv[1] + "/survey_%05i/paramfile_%05i.csv" % (i+1, i+1), big_survey = 0, **combs[i])

    if random.random() < 0.0:
        make_param_file(sys.argv[1] + "/survey_big_%05i/paramfile_%05i.csv" % (i+1, i+1), big_survey = 1, **combs[i])
        make_param_file(sys.argv[1] + "/survey_huge_%05i/paramfile_%05i.csv" % (i+1, i+1), big_survey = 2, **combs[i])



