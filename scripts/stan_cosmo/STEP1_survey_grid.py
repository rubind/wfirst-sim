import numpy as np
from subprocess import getoutput
import tqdm


def write_tier(f, total_survey_years, tier_name, tier_percent, exp_times, cadence, max_z):
    if tier_percent == 0:
        return 0

    cadence_steps = 365.24*2/cadence
    exp_times_w_overhead_seconds = sum(np.array(exp_times) + 52.)
    number_of_pointings_per_visit = (0.01*tier_percent*total_survey_years*86400*365.24)/(exp_times_w_overhead_seconds*cadence_steps)
    square_degrees = number_of_pointings_per_visit*0.281
    
    f.write("""__________________________________________,,,,,,
tier_name,%s,,,,,
tier_fraction_time,%.2f,,,,,
square_degrees,%.2f,,,,,
filters,R062,Z087,Y106,J129,H158,F184,
exp_times_per_dither,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,,
cadence,%i,%i,%i,%i,%i,%i,
dithers_per_filter,1,1,1,1,1,1,
trigger_redshifts,0,0.01,0.65,0.8,0.81,3,
trigger_fraction,1,0,0,0,0,0,
parallel_filters,H158,,,,,
max_SNe,1000000,
max_z,%.2f,\n""" % (tier_name, tier_percent/100., square_degrees,
                    exp_times[0], exp_times[1], exp_times[2],
                    exp_times[3], exp_times[4], exp_times[5],
                    cadence, cadence, cadence, cadence, cadence, cadence,
                    max_z))
    

def make_survey(total_survey_years, widepercent, medpercent, deeppercent):
    wd = "survey_grid/yr=%.3f_w=%03i_m=%03i_d=%03i" % (total_survey_years, widepercent, medpercent, deeppercent)
    getoutput("mkdir -p " + wd)

    f = open(wd + "/paramfile.csv", 'w')

    f.write("""total_survey_time,0.375,,,,,
maximum_trigger_fraction,0.9,,,,,
adjust_each_SN_exp_time,FALSE,,,,,
normalization_wavelength_range,5000,6000,,,,
shallow_SNR,13.5,,,,,
medium_SNR,23.2,,,,,
deep_SNR,38.6,,,,,
reference_SNR,45,,,,,
number_of_reference_dithers,8,,,,,
spectra_depth_and_phase,shallow  -10,medium -1,deep 1,shallow 5,shallow 10,
targeted_parallels,TRUE,,,,,
SN_rates,SN_rates.txt,,,,,
survey_duration,4,,,,,
grizY_30s_ground_depths,24.47,24.16,23.4,22.23,21.57,
hours_per_visit,30,,,,,
__________________________________________,,,,,,
slew_table,slew50.txt,,,,,
zodiacal_background,aldering.txt,,,,,
telescope_temperature,264,,,,,
PSFs,WebbPSF_WFC,,,,,
WFI_PSFs,WebbPSF_WFC,,,,,
interpixel_capacitance,0.02,,,,,
WFI_read_noise_floor,5,,,,,
WFI_read_noise_white,15,,,,,
WFI_dark_current,0.015,,,,,
WFI_pixel_scale,0.11,
IFU_read_noise_floor,4,,,,,
IFU_read_noise_white,15,,,,,
IFU_min_wave,4200,,,,,
IFU_max_wave,20000,,,,,
IFU_effective_area,IFU_effective_area_160720.txt,,,,,
IFU_resolution,IFU_R_160720.txt,,,,,
IFU_pixel_scale,0.05,,,,,
IFU_slice_in_pixels,3,,,,,
IFU_dark_current,0.003,,,,,
bad_pixel_rate,0.01,,,,,
__________________________________________,,,,,,
tier_name,Nearby,,,,,,
tier_fraction_time,0.,,,,,,
square_degrees,5000.,,,,,,
filters,g,r,i,z,
exp_times_per_dither,1,1,1,1,
cadence,4,4,4,4,
dithers_per_filter,1,1,1,1,
trigger_redshifts,0,0.01,0.65,0.8,0.81,3
trigger_fraction,1,0,0,0,0,0
parallel_filters,H158,,,,,
max_SNe,800,
max_z,0.1,
""")
    
    
    write_tier(f, total_survey_years = total_survey_years, tier_name = "Wide", tier_percent = widepercent, exp_times = [24.6, 31.4, 42.8, 61.8, 94, 175.4], cadence = 10, max_z = 1.0)
    write_tier(f, total_survey_years = total_survey_years, tier_name = "Medium", tier_percent = medpercent, exp_times = [152.9, 67.6, 75.3, 92.2, 187.9, 390.4], cadence = 5, max_z = 2.0)
    write_tier(f, total_survey_years = total_survey_years, tier_name = "Deep", tier_percent = deeppercent, exp_times = [152.9, 152.9, 235.4, 246.3, 336.7, 1017.6], cadence = 5, max_z = 2.5)

    f.close()

    pwd = getoutput("pwd")
    f = open(wd + "/run.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --partition=shared
#SBATCH --time=0-05:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile
export WFIRST=/home/drubin/wfirst-sim/
export WFIRST_SIM_DATA=/home/drubin/wfirst-sim-data/

""")
    f.write("cd " + pwd + "/" + wd + '\n')
    f.write("python ../../STEP1_simulate_survey.py paramfile.csv survey.pickle > log.txt\n")
    f.write("python ../../STEP2_Analytic_Fisher.py survey.pickle > fisher_log.txt\n")
    f.write("python ../../FoM.py comb_mat.fits > FoM.txt\n")
    f.close()

    print(getoutput("cd " + wd "\n sbatch run.sh"))
    

wide_exp_times = [] # RZYJHF, z=0.5
med_exp_times = [] # RZYJHF, z=1.0
deep_exp_times = [] # RZYJHF, z=1.7


getoutput("rm -fr survey_grid")

for widepercent in np.arange(0, 101, 5):
    for medpercent in np.arange(0, 101, 5):
        if widepercent + medpercent <= 100:
            deeppercent = 100 - (widepercent + medpercent)

            print("widepercent", widepercent, "medpercent", medpercent, "deeppercent", deeppercent)

            make_survey(total_survey_years = 0.375, widepercent = widepercent, medpercent = medpercent, deeppercent = deeppercent)
