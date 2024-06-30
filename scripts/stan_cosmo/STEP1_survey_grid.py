import numpy as np
from subprocess import getoutput
import tqdm
import sys


def write_tier(f, total_survey_years, tier_name, tier_percent, exp_times, filters, cadence, max_z):
    if tier_percent == 0:
        return 0

    cadence_steps = 365.24*2/cadence
    exp_times_w_overhead_seconds = sum(np.array(exp_times) + 52.)
    number_of_pointings_per_visit = (0.01*tier_percent*total_survey_years*86400*365.24)/(exp_times_w_overhead_seconds*cadence_steps)
    square_degrees = number_of_pointings_per_visit*0.281

    filter_names = [dict(R = "R062", Z = "Z087",
                         Y = "Y106", J = "J129",
                         H = "H158", F = "F184",
                         P = "P100")[item] for item in filters]
    
    f.write("""__________________________________________,,,,,,
tier_name,%s,,,,,
tier_fraction_time,%.4f,,,,,
square_degrees,%.2f,,,,,
filters,%s,
exp_times_per_dither,%s,
cadence,%s,
dithers_per_filter,%s,
max_SNe,1000000,
max_z,%.2f,\n""" % (tier_name, tier_percent/100., square_degrees,
                    ",".join(filter_names), ",".join(["%.2f" % item for item in exp_times]),
                    ",".join([str(cadence)]*len(filters)),
                    ",".join(["1"]*len(filters)),
                    max_z))
    

    
def make_survey(total_survey_years, widepercent, medpercent, deeppercent, widepercent_prism, deeppercent_prism, nnearby):
    wd = location + "/yr=%.3f_w=%03i_m=%03i_d=%03i_wp=%03i_dp=%03i_nnearby=%05i" % (total_survey_years, widepercent, medpercent, deeppercent, widepercent_prism, deeppercent_prism, nnearby)
    getoutput("mkdir -p " + wd)

    f = open(wd + "/paramfile.csv", 'w')

    square_degrees = 5000*(nnearby/800.)
    
    f.write("""total_survey_time,%.4f,,,,,
SN_rates,SN_rates.txt,,,,,
survey_duration,4,,,,,
total_FoV,0.82034*0.37628,
active_FoV,0.281*0.99,
grizY_30s_ground_depths,24.47,24.16,23.4,22.23,21.57,
hours_per_visit,30,,,,,
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
""" % (total_survey_years))
    
    if nnearby > 0:
        f.write("""__________________________________________,,,,,,
tier_name,Nearby,,,,,,
tier_fraction_time,0.,,,,,,
square_degrees,%i,,,,,,
filters,g,r,i,z,
exp_times_per_dither,1,1,1,1,
cadence,4,4,4,4,
dithers_per_filter,1,1,1,1,
trigger_redshifts,0,0.01,0.65,0.8,0.81,3
trigger_fraction,1,0,0,0,0,0
parallel_filters,H158,,,,,
max_SNe,%i,
max_z,0.1,
""" % (square_degrees, nnearby))



        
    write_tier(f, total_survey_years = total_survey_years, tier_name = "WideNoPrism", tier_percent = widepercent*(1. - widepercent_prism/100.), exp_times = [24.6, 31.4, 42.8, 61.8, 94, 175.4], filters = "RZYJHF", cadence = 10, max_z = 1.0)
    write_tier(f, total_survey_years = total_survey_years, tier_name = "WidePrism", tier_percent = widepercent*widepercent_prism/100., exp_times = [24.6, 31.4, 42.8, 61.8, 94, 175.4, 1800.], filters = "RZYJHFP", cadence = 10, max_z = 1.0)
    write_tier(f, total_survey_years = total_survey_years, tier_name = "MediumNoPrism", tier_percent = medpercent, exp_times = [152.9, 67.6, 75.3, 92.2, 187.9, 390.4], cadence = 5, filters = "RZYJHF", max_z = 2.0)
    write_tier(f, total_survey_years = total_survey_years, tier_name = "DeepNoPrism", tier_percent = deeppercent*(1. - deeppercent_prism/100.), exp_times = [152.9, 152.9, 235.4, 246.3, 336.7, 1017.6], cadence = 5, filters = "RZYJHF", max_z = 2.5)
    write_tier(f, total_survey_years = total_survey_years, tier_name = "DeepPrism", tier_percent = deeppercent*deeppercent_prism/100., exp_times = [152.9, 152.9, 235.4, 246.3, 336.7, 1017.6, 3600.], cadence = 5, filters = "RZYJHFP", max_z = 2.5)

    f.close()

    memory_needed = 12 + 12*(deeppercent_prism + widepercent_prism > 0)
    
    pwd = getoutput("pwd")
    f = open(wd + "/run.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --partition=shared
#SBATCH --time=0-05:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=%iG # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile
export WFIRST=/home/drubin/wfirst-sim/
export WFIRST_SIM_DATA=/home/drubin/wfirst-sim-data/
pip install sncosmo
pip install sep

""" % memory_needed)
    f.write("cd " + pwd + "/" + wd + '\n')
    f.write("python $WFIRST/scripts/stan_cosmo/STEP1_simple_survey.py paramfile.csv survey.pickle > log.txt\n")
    f.write("python $WFIRST/scripts/stan_cosmo/STEP2_Analytic_Fisher.py survey.pickle > fisher_log.txt\n")
    f.write("python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat.fits > FoM.txt\n")
    f.close()

    print(getoutput("cd " + wd + "\n sbatch run.sh"))
    

wide_exp_times = [] # RZYJHF, z=0.5
med_exp_times = [] # RZYJHF, z=1.0
deep_exp_times = [] # RZYJHF, z=1.7

grid_type = sys.argv[1]
location = sys.argv[2]

getoutput("rm -fr " + location)

if grid_type == "tier_fraction":
    for widepercent in np.arange(0, 101, 5):
        for medpercent in np.arange(0, 101, 5):
            if widepercent + medpercent <= 100:
                deeppercent = 100 - (widepercent + medpercent)
                
                print("widepercent", widepercent, "medpercent", medpercent, "deeppercent", deeppercent)
                
                make_survey(total_survey_years = 0.5, widepercent = widepercent, medpercent = medpercent, deeppercent = deeppercent, nnearby = 800, widepercent_prism = 0, deeppercent_prism = 25)

elif grid_type == "total_time":
    for total_survey_years in np.arange(0.05, 1.01, 0.025):
        make_survey(total_survey_years = total_survey_years, widepercent = 70, medpercent = 0, deeppercent = 30, nnearby = 800)

elif grid_type == "nnearby":
    for nnearby in np.arange(0, 3001, 200):
        make_survey(total_survey_years = 0.375, widepercent = 70, medpercent = 0, deeppercent = 30, nnearby = int(nnearby))

elif grid_type == "prism_fraction":
    for widepercent_prism in np.arange(0, 21, 1):
        for deeppercent_prism in np.arange(0, 101, 5):
            deeppercent = 30
            widepercent = 70

            print("widepercent_prism", widepercent_prism, "deeppercent_prism", deeppercent_prism, "widepercent", widepercent, "deeppercent", deeppercent)
            make_survey(total_survey_years = 0.5, widepercent = widepercent, medpercent = 0, deeppercent = deeppercent, nnearby = 800, widepercent_prism = widepercent_prism, deeppercent_prism = deeppercent_prism)
            


else:
    print("Unknown grid type! want: tier_fraction total_time prism_fraction or nnearby")
