import numpy as np
from subprocess import getoutput
import tqdm
import sys
from FileRead import readcol, format_file
import os
WFIRST = os.environ["WFIRST"]

minimum_exp_time = 60.
minimum_rest_wave_target = 0.37 # microns


def get_per_visit_exposure_time(band, redshift, cadence):
    f = open(WFIRST + "/scripts/pixel-level/StoN_vs_exptime_%s_TTel_264.0_WebbPSF_WFC_aldering_dark=0.015_rnf=5.0_rnw=16.0_SNRtarg=10.00@1.00_addhost=1_mod=SALT3.NIR_WAVEEXT_zlots=1.0.txt" % (band), 'r')
    lines = f.read().split('\n')
    f.close()

    rest_wave_micron = dict(R062 = 0.62, Z087 = 0.87, Y106 = 1.06, J129 = 1.29, H158 = 1.58, F184 = 1.84, K213 = 2.13)[band]

    redshift = min(redshift, rest_wave_micron/minimum_rest_wave_target - 1.)
    redshift = np.around(redshift, 1)

    for line in lines:
        if line.count("50.0: ") == 1:
            if line.count("redshift=%.1f" % redshift) == 1:
                if line.count("key=%.2f" % (10.*np.sqrt(   cadence/(2.5*(1. + redshift))   ))) == 1:
                    exp_in_five_days = float(line.split("exp=")[-1])
                    exp_per_visit = (exp_in_five_days/5.)*cadence
                    exp_per_visit = max(exp_per_visit, minimum_exp_time)
                    return exp_per_visit
                
    assert 0, "Error reading " + str(band) + " " + str(redshift) + " " + str(cadence)


f = open("all_exp_times.txt", 'w')
tmp_cadences = np.arange(2, 15.1, 1)

for band in ["R062", "Z087", "Y106", "J129", "H158", "F184"]:
    f.write("Band Redshift " + " ".join(["cad=%.1f" % item for item in tmp_cadences]) + '\n')
    for redshift in np.arange(0.2, 2.51, 0.1):
        f.write(band + " z=%.2f" % redshift)
        for cadence in tmp_cadences:
            f.write(" %.2f" % get_per_visit_exposure_time(band = band, redshift = redshift, cadence = cadence))
        f.write('\n')
f.close()
format_file("all_exp_times.txt")

f = open("all_exp_times.txt", 'r')
lines = f.read()
f.close()

lines = "#Notes:\n" + ("#No exposure times under %.1fs are allowed.\n" % minimum_exp_time) + ("#The redshift targeted cannot go above a redshift where the filter goes bluer than rest-frame %.3f microns, e.g., z=%.2f for R062.\n" % (minimum_rest_wave_target, 0.62/minimum_rest_wave_target - 1.)) + "#In the Poisson regime, exposure time scales linearly with cadence when targeting a specific redshift, as it should.\n" + lines

f = open("all_exp_times.txt", 'w')
f.write(lines)
f.close()





def get_exp_time_dict(redshift, cadence):
    exp_time_dict = {}
    
    for band in ["R062", "Z087", "Y106", "J129", "H158", "F184"]:
        exp_time_dict[band[0]] = get_per_visit_exposure_time(band = band, redshift = redshift, cadence = cadence)
    return exp_time_dict
    
def get_square_degrees(cadence, filters, exp_times, tier_percent, total_survey_years, prism_cadence = 5.):
    if filters.count("g"):
        survey_duration = 0.6
    else:
        assert filters.count("r") == 0
        assert filters.count("i") == 0
        assert filters.count("z") == 0

        survey_duration = 2.0

    cadence_steps = 365.24*survey_duration/cadence
    

    exp_times_roman_w_overhead = []
    for i in range(len(filters)):
        if "grizy".count(filters[i]) == 0:
            if filters[i].count("P") == 0:
                exp_times_roman_w_overhead.append(exp_times[i] + 53.)
            else:
                exp_times_roman_w_overhead.append((exp_times[i] + 53.)*cadence/prism_cadence)


    print("exp_times_roman_w_overhead", exp_times_roman_w_overhead)

    total_w_overhead_seconds = sum(exp_times_roman_w_overhead)
    number_of_pointings_per_visit = (0.01*tier_percent*total_survey_years*86400*365.24)/(total_w_overhead_seconds*cadence_steps)
    square_degrees = number_of_pointings_per_visit*0.82034*0.37628

    return square_degrees, survey_duration


def write_tier(f, total_survey_years, tier_name, tier_percent, exp_times, filters, cadence, max_z, prism_cadence = 5.):
    if tier_percent == 0:
        return 0

    assert len(filters) == len(exp_times)
        
    square_degrees, survey_duration = get_square_degrees(cadence = cadence, filters = filters, exp_times = exp_times, tier_percent = tier_percent, total_survey_years = total_survey_years, prism_cadence = prism_cadence)

    
    filter_names = [dict(R = "R062", Z = "Z087",
                         Y = "Y106", J = "J129",
                         H = "H158", F = "F184",
                         P = "P100", u = "u", g = "g", r = "r", i = "i", z = "z", y = "y")[item] for item in filters]


    cadences = []
    cadence_offsets = []
    
    for i, filt in enumerate(filters):
        if "ugrizy".count(filt):
            cadences.append("4")
            cadence_offsets.append("0")
        elif filt == "P":
            cadences.append(str(prism_cadence))
            cadence_offsets.append("0")
        else:
            cadences.append(str(cadence))
            cadence_offsets.append(str(np.around((i%2)*cadence/2.)))
            
    
    f.write("""__________________________________________,,,,,,
tier_name,%s,,,,,
tier_fraction_time,%.4f,,,,,
square_degrees,%.2f,,,,,
filters,%s,
exp_times_per_dither,%s,
cadences,%s,
cadence_offsets,%s,
dithers_per_filter,%s,
max_SNe,1000000,
max_z,%.2f,\n""" % (tier_name, tier_percent/100., square_degrees,
                    ",".join(filter_names), ",".join(["%.2f" % item for item in exp_times]),
                    ",".join(cadences),
                    ",".join(cadence_offsets),
                    ",".join(["1"]*len(filters)),
                    max_z))
    

    
def make_survey(total_survey_years, widepercent_imaging, medpercent_imaging, widepercent_prism, medpercent_prism, deeppercent_prism,
                nnearby, wide_filts, med_filts, deep_filts, SN_number_poisson,
                #exp_times_dict_wide,
                #exp_times_dict_med,
                #exp_times_dict_deep,
                suffix = "0", wd = "",
                add_Rubin_only_tier = 0,
                wide_cadence = 10, med_cadence = 5, deep_cadence = 5,
                wide_ztarg = 0.5, med_ztarg = 1.0, deep_ztarg = 1.7,
                SN_rates = "SN_rate1.txt", SNRMax = 0, wide_rubin = 0, model_res_list = [9]):

    
    deeppercent_imaging = 100 - (widepercent_imaging + medpercent_imaging + widepercent_prism + medpercent_prism + deeppercent_prism)

    if widepercent_imaging < 0 or widepercent_imaging > 100:
        return 0

    if medpercent_imaging < 0 or medpercent_imaging > 100:
        return 0
    
    if deeppercent_imaging < 0 or deeppercent_imaging > 100:
        return 0

    

    if wide_rubin:
        wide_filts += "griz"

    if len(wide_filts) < 3:
        return 0

    if len(med_filts) < 3:
        return 0

    if len(deep_filts) < 3:
        return 0

    if med_ztarg <= wide_ztarg + 0.001:
        return 0
    if deep_ztarg <= wide_ztarg + 0.001:
        return 0
    if deep_ztarg <= med_ztarg + 0.001:
        return 0



    exp_times_dict_wide = get_exp_time_dict(redshift = wide_ztarg, cadence = wide_cadence)
    exp_times_dict_med = get_exp_time_dict(redshift = med_ztarg, cadence = med_cadence)
    exp_times_dict_deep = get_exp_time_dict(redshift = deep_ztarg, cadence = deep_cadence)

        
    if not "P" in exp_times_dict_wide:
        exp_times_dict_wide["P"] = 900
    if not "P" in exp_times_dict_med:
        exp_times_dict_med["P"] = 1800
    if not "P" in exp_times_dict_deep:
        exp_times_dict_deep["P"] = 3600


    exp_times = [exp_times_dict_wide[item] for item in wide_filts]
    sq_wide_imaging, NA = get_square_degrees(cadence = wide_cadence, filters = wide_filts, exp_times = exp_times, tier_percent = widepercent_imaging, total_survey_years = total_survey_years)
    sq_wide_prism, NA = get_square_degrees(cadence = wide_cadence, filters = "P", exp_times = [exp_times_dict_wide["P"]], tier_percent = widepercent_prism, total_survey_years = total_survey_years)

    exp_times = [exp_times_dict_med[item] for item in med_filts]
    sq_med_imaging, NA = get_square_degrees(cadence = med_cadence, filters = med_filts, exp_times = exp_times, tier_percent = medpercent_imaging, total_survey_years = total_survey_years)
    sq_med_prism, NA = get_square_degrees(cadence = med_cadence, filters = "P", exp_times = [exp_times_dict_med["P"]], tier_percent = medpercent_prism, total_survey_years = total_survey_years)

    exp_times = [exp_times_dict_deep[item] for item in deep_filts]
    sq_deep_imaging, NA = get_square_degrees(cadence = deep_cadence, filters = deep_filts, exp_times = exp_times, tier_percent = deeppercent_imaging, total_survey_years = total_survey_years)
    sq_deep_prism, NA = get_square_degrees(cadence = deep_cadence, filters = "P", exp_times = [exp_times_dict_deep["P"]], tier_percent = deeppercent_prism, total_survey_years = total_survey_years)
    
    
    if wd == "":
        wd = location + "/yr=%.3f_wi=%03i_mi=%03i_di=%03i_wp=%03i_mp=%03i_dp=%03i_nnearby=%05i_%s+%s+%s_cad=%02i+%02i+%02i_PN=%i_Ronly=%i_ztarg=%.1f+%.1f+%.1f_sq=%.2g+%.2g+%.2g_%s" % (total_survey_years, widepercent_imaging, medpercent_imaging, deeppercent_imaging,
                                                                                                                                                                      widepercent_prism, medpercent_prism, deeppercent_prism,
                                                                                                                                                                      nnearby, wide_filts, med_filts, deep_filts, wide_cadence, med_cadence, deep_cadence,
                                                                                                                                                                      SN_number_poisson, add_Rubin_only_tier,
                                                                                                                                                                      wide_ztarg, med_ztarg, deep_ztarg,
																				      sq_wide_imaging, sq_med_imaging, sq_deep_imaging,
                                                                                                                                                                      suffix)
    else:
        wd = location + "/" + wd

        
    getoutput("mkdir -p " + wd)

    f = open(wd + "/paramfile.csv", 'w')

    square_degrees = 5000*(nnearby/800.)
    
    f.write("""total_survey_time,%.4f,,,,,
SN_rates,%s,,,,,
survey_duration,2.5,,,,,
total_FoV,0.82034*0.37628,
active_FoV,0.281*0.99,
grizY_30s_ground_depths,24.47,24.16,23.4,22.23,21.57,
hours_per_visit,30,,,,,
sn_model,SALT3.NIR_WAVEEXT,
SN_number_poisson,%i,
__________________________________________,,,,,,
slew_table,slew50.txt,,,,,
zodiacal_background,aldering.txt,,,,,
telescope_temperature,264,,,,,
PSFs,Prism_3Element_PSF_SCA08,,,,,
WFI_PSFs,WebbPSF_WFC,,,,,
interpixel_capacitance,0.02,,,,,
WFI_read_noise_floor,5,,,,,
WFI_read_noise_white,16,,,,,
WFI_dark_current,0.015,,,,,
WFI_pixel_scale,0.11,
IFU_read_noise_floor,5,,,,,
IFU_read_noise_white,16,,,,,
IFU_min_wave,7500,,,,,
IFU_max_wave,18000,,,,,
IFU_effective_area,P100.txt,,,,,
IFU_resolution,Prism_R.txt,,,,,
IFU_pixel_scale,0.11,,,,,
IFU_slice_in_pixels,5,,,,,
IFU_dark_current,1.2,,,,,
bad_pixel_rate,0.01,,,,,
""" % (total_survey_years, SN_rates, SN_number_poisson))
    
    if nnearby > 0:
        f.write("""__________________________________________,,,,,,
tier_name,Nearby,,,,,,
tier_fraction_time,0.,,,,,,
square_degrees,%i,,,,,,
filters,g,r,i,z,
exp_times_per_dither,1,1,1,1,
cadences,4,4,4,4,
cadence_offsets,0,0,0,0,
dithers_per_filter,1,1,1,1,
max_SNe,%i,
max_z,0.1,
""" % (square_degrees, nnearby))

    if add_Rubin_only_tier:
        f.write("""__________________________________________,,,,,,
tier_name,RubinDDF,,,,,,
tier_fraction_time,0.,,,,,,
square_degrees,9.6*8,,,,,,
filters,g,r,i,z,
exp_times_per_dither,30,30,30,30,
cadences,4,4,4,4,
cadence_offsets,0,0,0,0,
dithers_per_filter,1,1,1,1,
max_SNe,100000,
max_z,0.5,
""")


    imaging_bins = [0.,
                    sq_deep_imaging,
                    sq_deep_imaging + sq_med_imaging,
                    sq_wide_imaging + sq_med_imaging + sq_deep_imaging]

    prism_bins = [0.,
                  sq_deep_prism,
                  sq_deep_prism + sq_med_prism,
                  sq_wide_prism + sq_med_prism + sq_deep_prism,
                  1e10]

    print("imaging_bins", imaging_bins)
    print("prism_bins", prism_bins)

    
    if prism_bins[-2] > imaging_bins[-1]:
        print("Prism larger than imaging!")
        return 0


    all_percent = 0.
    for ind_im, im_name in enumerate(["DeepImaging", "MediumImaging", "WideImaging"]):
        for ind_pr, pr_name in enumerate(["DeepPrism", "MediumPrism", "WidePrism", "NoPrism"]):
            tier_min = max(imaging_bins[ind_im], prism_bins[ind_pr])
            tier_max = min(imaging_bins[ind_im + 1], prism_bins[ind_pr + 1])

            tier_sq = tier_max - tier_min
            if tier_sq > 0:
                overlap_imaging = tier_sq/([sq_deep_imaging, sq_med_imaging, sq_wide_imaging, 1e10][ind_im])
                overlap_prism =	tier_sq/([sq_deep_prism, sq_med_prism, sq_wide_prism, 1e10][ind_pr])

                
                tier_percent = [deeppercent_imaging, medpercent_imaging, widepercent_imaging][ind_im]*overlap_imaging + [deeppercent_prism, medpercent_prism, widepercent_prism, 0.][ind_pr]*overlap_prism
                all_percent += tier_percent

                filters = [deep_filts, med_filts, wide_filts][ind_im]
                exp_times = [[exp_times_dict_deep, exp_times_dict_med, exp_times_dict_wide][ind_im][item] for item in filters]

                if pr_name != "NoPrism":
                    filters += "P"
                    exp_times.append([exp_times_dict_deep, exp_times_dict_med, exp_times_dict_wide][ind_pr]["P"])

                
                write_tier(f, total_survey_years = total_survey_years,
                           tier_name = im_name + pr_name, tier_percent = tier_percent, exp_times = exp_times, filters = filters,
                           cadence = [deep_cadence, med_cadence, wide_cadence][ind_im],
                           prism_cadence = 5.,
                           max_z = [3.0, 3.0, 3.0][ind_im])
    f.close()
    assert np.isclose(all_percent, 100)
            

    memory_needed = 16 + 16*(deeppercent_prism + medpercent_prism + widepercent_prism > 0)

    
    pwd = getoutput("pwd")
    f = open(wd + "/run.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --partition=shared
#SBATCH --time=0-08:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=""" + str(memory_needed) + """G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile
export WFIRST=/home/drubin/wfirst-sim/
export WFIRST_SIM_DATA=/home/drubin/wfirst-sim-data/
pip install sncosmo
pip install sep

""")
    f.write("cd " + pwd + "/" + wd + '\n')
    f.write("python $WFIRST/scripts/stan_cosmo/STEP1_simple_survey.py paramfile.csv survey.pickle > log.txt\n")

    for model_res in model_res_list:
        f.write("python $WFIRST/scripts/stan_cosmo/STEP2_Analytic_Fisher.py survey.pickle --SNRMax %i --model_res %i  --train 1 --calib 1 --color_scatter_opt 0.04 --color_scatter_nir 0.02 --gray_disp 0.08 > fisher_log.txt\n" % (SNRMax, model_res))
        f.write("python $WFIRST/scripts/stan_cosmo/STEP1A_plot_survey.py survey.pickle\n")
        f.write("python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat.fits > FoM_model_res=%02i.txt\n" % model_res)
        f.write("python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat_no_model.fits > FoM_no_model_%02i.txt\n" % model_res)
        f.write("python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat_stat_only.fits > FoM_stat_only_%02i.txt\n" % model_res)
    f.close()

    print(getoutput("cd " + wd + "\n sbatch run.sh"))
    return 1


def make_survey_wrap(*args, **kwargs):
    for realization in range(n_real):
        kwargs["suffix"] = "%02i" % realization
        make_survey(*args, **kwargs)


def read_csv():
    # 2TIER_PRISM25_4bands a00-t00-z00,SHALLOW,RZYJ,0.5,91,25.72,242,3,R:30.0;Z:30.0;Y:33.0;J:47.2,DEEP,YJHF,1.2,22,6.43,146,5,Y:101.3;J:109.1;H:201.0;F:471.5

    [survey_name, tier_name, wide_filts, NA, NA, wide_square_degrees, wide_nvisit, wide_cadence, wide_exptimes,
     tier_name, deep_filts, NA, NA, deep_square_degrees, deep_nvisit, deep_cadence, deep_exptimes] = readcol("params.csv", "aaa,ff,fff,a,,aa,ff,fff,a", splitchar = ",")

    for i in range(len(survey_name)):
        exp_times_dict_wide = {}
        exp_times_dict_deep = {}
        
        for item in wide_exptimes[i].split(";"):
            exp_times_dict_wide[item.split(":")[0]] = float(item.split(":")[1])

        for item in deep_exptimes[i].split(";"):
            exp_times_dict_deep[item.split(":")[0]] = float(item.split(":")[1])

        print(exp_times_dict_wide)
        print(exp_times_dict_deep)

        tot_wide = sum([exp_times_dict_wide[item] + 55. for item in exp_times_dict_wide])
        tot_deep = sum([exp_times_dict_deep[item] + 55. for item in exp_times_dict_deep])

        tot_wide *= wide_nvisit[i]*wide_square_degrees[i]/0.28
        tot_deep *= deep_nvisit[i]*deep_square_degrees[i]/0.28

        print("total ", (tot_wide + tot_deep)/(86400*365.24))
        
        make_survey(total_survey_years = 0.375, widepercent_imaging = 100.*tot_wide/(tot_wide + tot_deep), medpercent = 0., widepercent_prism = 0., medpercent_prism = 0., deeppercent_prism = 0., nnearby = 800.,
                    wide_filts = wide_filts[i], med_filts = "", deep_filts = deep_filts[i], exp_times_dict_wide = exp_times_dict_wide, exp_times_dict_deep = exp_times_dict_deep,
                    wd = survey_name[i].replace(" ", "_"), wide_cadence = wide_cadence[i], deep_cadence = deep_cadence[i], SN_rates = "SN_rates_powerlaw.txt", SNRMax = 1)
        
        
wide_exp_times = [] # RZYJHF, z=0.5
med_exp_times = [] # RZYJHF, z=1.0
deep_exp_times = [] # RZYJHF, z=1.7

percent_prism10000 = []
for i in range(10000):
    this_percent = 10000
    while this_percent >= 99.9:
        this_percent = np.random.exponential()*30
    percent_prism10000.append(this_percent)
    

grid_vals = dict(widepercent_imaging = np.arange(0, 201, 5),
                 medpercent_imaging = np.arange(0, 101, 5),
                 deeppercent_imaging = np.arange(0, 201, 5),
                 total_survey_years = [0.5],
                 nnearby = [800],
                 add_Rubin_only_tier = [1],
                 wide_filts = ["RZYJHF"]*3 + ["RZYJH"]*3 + ["RZYJ", "RZJH"]*3 + ["RZY", "RZJ", "RZH", "ZYJ", "ZYJH", "ZJH", "ZHF", "ZYH"],
                 med_filts = ["RZYJHF", "ZYJHF", "YJHF", "RZYJH", "RZYJ", "ZYJH"],
	         deep_filts = ["RZYJHF", "ZYJHF", "YJHF", "RZYJH", "ZYJH"],
                 widepercent_prism = np.arange(0, 41, 2), #[0],#np.arange(0, 41, 2),
                 medpercent_prism = np.arange(0, 41, 2), #[0],#np.arange(0, 41, 2),
                 deeppercent_prism = np.arange(0, 41, 2), #[0],#np.arange(0, 41, 2),
                 percent_prism = np.array(percent_prism10000), #np.arange(0, 16, 1.), #np.array(percent_prism10000),
                 SN_number_poisson = [0])
                 

grid_vals = dict(widepercent_imaging = np.arange(0, 101, 5),
                 medpercent_imaging = np.arange(0, 101, 5),
                 deeppercent_imaging = np.arange(0, 101, 5),
                 total_survey_years = [0.5],
                 nnearby = [800],
                 add_Rubin_only_tier = [1],
                 wide_filts = ["RZJRHY"]*4 + ["RZYJHF"]*3 + ["RZYJH"]*3 +  ["RZYJ", "RZJH"]*2 + ["RZY", "RZJ", "RZH", "ZYJ", "ZYJH", "ZJH", "ZHF", "ZYH"],
                 med_filts = ["RZYJHF", "ZYJHF", "YJHF", "RZYJH", "RZYJ", "ZYJH"],
                 deep_filts = ["RZYJHF", "ZYJHF"],
                 wide_cadence = [10],
                 wide_ztarg = np.arange(0.3, 1.1, 0.1), med_ztarg = np.arange(0.8, 1.3, 0.1), deep_ztarg = np.arange(1.0, 2.4, 0.1),
                 med_cadence = [4, 5, 6, 7, 8, 9, 10, 12, 15],
                 deep_cadence = [4, 5, 6, 7, 8, 9, 10, 12, 15],
                 widepercent_prism = np.arange(0, 41, 2), #[0],#np.arange(0, 41, 2),
                 medpercent_prism = np.arange(0, 41, 2), #[0],#np.arange(0, 41, 2),
                 deeppercent_prism = np.arange(0, 41, 2), #[0],#np.arange(0, 41, 2),
                 percent_prism = np.array(percent_prism10000), #np.arange(0, 16, 1.), #np.array(percent_prism10000),                                                                                                                                 
                 SN_number_poisson = [0])



#exp_times_dict_wide = dict(R = 24.6, Z = 31.4, Y = 42.8, J = 61.8, H = 94, F = 175.4, P = 900., g = 30., r = 30., i = 30., z = 30.) # Note, cadence is 10 for imaging

#exp_times_dict_wide = dict(R = 60., Z = 60., Y = 60., J = 61.8, H = 94, F = 175.4, P = 900., g = 30., r = 30., i = 30., z = 30.) # Minimum 60s
#exp_times_dict_med = dict(R = 152.9, Z = 67.6, Y = 75.3, J = 92.2, H = 187.9, F = 390.4, P = 1800.)
#exp_times_dict_deep = dict(R = 152.9, Z = 152.9, Y = 235.4, J = 246.3, H = 336.7, F = 1017.6, P = 3600.)



grid_type = sys.argv[1]
location = sys.argv[2]
n_real = int(sys.argv[3])

getoutput("rm -fr " + location)




if grid_type == "tier_fraction":
    for widepercent_imaging in grid_vals["widepercent_imaging"]:
        for medpercent_imaging in grid_vals["medpercent_imaging"]:
            deeppercent_imaging = 100 - (widepercent_imaging + medpercent_imaging)
            
            print("widepercent_imaging", widepercent_imaging, "medpercent_imaging", medpercent_imaging, "deeppercent", deeppercent_imaging)

            for add_Rubin_only_tier in [0, 1]:
                make_survey(total_survey_years = 0.5, widepercent_imaging = widepercent_imaging, medpercent_imaging = medpercent_imaging, nnearby = 800, widepercent_prism = 0, medpercent_prism = 10, deeppercent_prism = 10,
                            wide_filts = "RZYJ", med_filts = "ZYJH", deep_filts = "ZYJHF", SN_number_poisson = 0, add_Rubin_only_tier = add_Rubin_only_tier,
                            exp_times_dict_wide = exp_times_dict_wide, exp_times_dict_med = exp_times_dict_med, exp_times_dict_deep = exp_times_dict_deep)

elif grid_type == "poisson":
    for SN_number_poisson in [0, 1]:
        make_survey_wrap(total_survey_years = 0.5, widepercent = 60, medpercent = 10, nnearby = 800, widepercent_prism = 0, medpercent_prism = 0, deeppercent_prism = 25,
                         wide_filts = "RZYJ", med_filts = "ZYJH", deep_filts = "ZYJHF", SN_number_poisson = SN_number_poisson)

elif grid_type == "total_time":
    for total_survey_years in np.arange(0.05, 1.01, 0.025):
        make_survey(total_survey_years = total_survey_years, widepercent = 70, medpercent = 0, deeppercent = 30, nnearby = 800)

elif grid_type == "nnearby":
    for nnearby in np.arange(0, 3001, 200):
        make_survey(total_survey_years = 0.375, widepercent = 70, medpercent = 0, deeppercent = 30, nnearby = int(nnearby))

elif grid_type == "filt_choice":
    for wide_filts in ["RZYJHF", "RZYJH", "RZYJ", "RZY", "RZJ", "RZH"]:
        for deep_filts in ["RZYJHF", "ZYJHF", "YJHF", "RZYJH", "ZYJH"]:
            make_survey_wrap(total_survey_years = 0.5, widepercent = 70, medpercent = 0, deeppercent = 30, nnearby = 800, widepercent_prism = 0, deeppercent_prism = 25,
                             wide_filts = wide_filts, med_filts = "RZYJHF", deep_filts = deep_filts)



elif grid_type == "prism_fraction":
    for widepercent_prism in np.arange(0, 21, 1):
        for deeppercent_prism in np.arange(0, 101, 5):
            deeppercent = 30
            widepercent = 70

            print("widepercent_prism", widepercent_prism, "deeppercent_prism", deeppercent_prism, "widepercent", widepercent, "deeppercent", deeppercent)
            make_survey(total_survey_years = 0.5, widepercent = widepercent, medpercent = 0, deeppercent = deeppercent, nnearby = 800, widepercent_prism = widepercent_prism, deeppercent_prism = deeppercent_prism)


elif grid_type == "prism_exp":
    for widepercent_imaging in [50.]:#[40., 50., 60.]:
        for widepercent_prism in tqdm.tqdm(np.arange(0., 25., 1.)):
            for deeppercent_prism in np.arange(0., 25., 1.):
                
                                
                make_survey_wrap(total_survey_years = 0.5, widepercent_imaging = widepercent_imaging, medpercent_imaging = 0., nnearby = 800, widepercent_prism = widepercent_prism, medpercent_prism = 0., deeppercent_prism = deeppercent_prism,
                                 wide_filts = "RZYJ", med_filts = "ZYJHF", deep_filts = "ZYJHF",
                                 exp_times_dict_wide = exp_times_dict_wide, exp_times_dict_med = exp_times_dict_med, exp_times_dict_deep = exp_times_dict_deep,
                                 SN_number_poisson = 0)

                
elif grid_type == "model_res":
    make_survey_wrap(total_survey_years = 0.375, widepercent = 50., medpercent = 0.0, nnearby = 800, widepercent_prism = 0, medpercent_prism = 0.,
                     deeppercent_prism = 0.0,
                     wide_filts = "RZYJ", med_filts = "ZYJHF", deep_filts = "ZYJHF",
                     exp_times_dict_wide = exp_times_dict_wide, exp_times_dict_med = exp_times_dict_med, exp_times_dict_deep = exp_times_dict_deep,
                     SN_number_poisson = 0, model_res_list = np.arange(7, 19))

elif grid_type == "cadence":
    for wide_cadence in np.arange(1., 16., 1.):
        for deep_cadence in np.arange(1., 16., 1.):
            exp_times_dict_wide = dict(R = get_per_visit_exposure_time(band = "R062", redshift = 0.5, cadence = wide_cadence),
                                       Z = get_per_visit_exposure_time(band = "Z087", redshift = 0.5, cadence = wide_cadence),
                                       Y = get_per_visit_exposure_time(band = "Y106", redshift = 0.5, cadence = wide_cadence),
                                       J = get_per_visit_exposure_time(band = "J129", redshift = 0.5, cadence = wide_cadence))
                                       
            exp_times_dict_deep = dict(Y = get_per_visit_exposure_time(band = "Y106", redshift = 1.7, cadence = deep_cadence),
                                       J = get_per_visit_exposure_time(band = "J129", redshift = 1.7, cadence = deep_cadence),
                                       H = get_per_visit_exposure_time(band = "H158", redshift = 1.7, cadence = deep_cadence),
                                       F = get_per_visit_exposure_time(band = "F184", redshift = 1.7, cadence = deep_cadence))

            
            make_survey_wrap(total_survey_years = 0.5, widepercent_imaging = 50., medpercent_imaging = 0.0, nnearby = 800, widepercent_prism = 0, medpercent_prism = 0.,
                             deeppercent_prism = 0.0,
                             wide_cadence = wide_cadence, med_cadence = 5, deep_cadence = deep_cadence,
                             wide_filts = "RZYJ", med_filts = "ZYJHF", deep_filts = "YJHF",
                             exp_times_dict_wide = exp_times_dict_wide, exp_times_dict_med = exp_times_dict_med, exp_times_dict_deep = exp_times_dict_deep,
                             SN_number_poisson = 0)
    

elif grid_type == "random":
    good_surveys = 0
    while good_surveys < n_real:
        these_pars = {}
        for key in grid_vals:
            these_pars[key] = np.random.choice(grid_vals[key])

        
        #tot_norm = these_pars["widepercent_imaging"] + these_pars["medpercent_imaging"] + these_pars["deeppercent_imaging"] + these_pars["widepercent_prism"] + these_pars["medpercent_prism"] + these_pars["deeppercent_prism"]

        im_norm = (100 - these_pars["percent_prism"])/(these_pars["widepercent_imaging"] + these_pars["medpercent_imaging"] + these_pars["deeppercent_imaging"])
        pr_norm = these_pars["percent_prism"]/(these_pars["widepercent_prism"] + these_pars["medpercent_prism"] + these_pars["deeppercent_prism"])
        
        
        for key1 in ["wide", "med", "deep"]:
            these_pars[key1 + "percent_imaging"] *= im_norm
            these_pars[key1 + "percent_imaging"] = int(np.around(these_pars[key1 + "percent_imaging"]))

            these_pars[key1 + "percent_prism"] *= pr_norm
            these_pars[key1 + "percent_prism"] = int(np.around(these_pars[key1 + "percent_prism"]))

        del these_pars["deeppercent_imaging"]
        del these_pars["percent_prism"]
        
        for wide_rubin in [0]: #,1]:
            these_pars["wide_rubin"] = wide_rubin
            print(these_pars)
            good_surveys += make_survey(**these_pars)
            print("good_surveys", good_surveys)

elif grid_type == "read_csv":
    read_csv()
    
        
else:
    print("Unknown grid type! want: tier_fraction total_time prism_fraction filt_choice nnearby read_csv poisson prism_exp model_res cadence or random")
