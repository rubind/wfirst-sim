from numpy import *
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import commands
import sys

def get_exptime(filt, max_z):
    if max_z == 0.8:
        return dict(Z087 = 50.7, Y106 = 55.0, J129 = 86.8, H158 = 137.6, F184 = 323.4)[filt]
    if max_z == 1.2:
        return dict(Z087 = 108.4, Y106 = 111.4, J129 = 132.4, H158 = 239.0, F184 = 639.5)[filt]
    if max_z == 1.7:
        return dict(Y106 = 243.3, J129 = 268.8, H158 = 361.6, F184 = 1003.7)[filt]
    

def get_filters(include_f184_for_08, max_z, spec_SNR):
    if spec_SNR == None:
        parfilters = []

    if max_z == 0.8:
        if spec_SNR == None:
            filters = ["Z087", "Y106", "J129", "H158"] + ["F184"]*include_f184_for_08
        else:
            filters = ["Z087", "Y106", "J129"]
            parfilters = ["H158"] + ["F184"]*include_f184_for_08

    elif max_z == 1.2:
        if spec_SNR == None:
            filters = ["Z087", "Y106", "J129", "H158", "F184"]
        else:
            filters = ["Z087", "Y106", "J129"]
            parfilters = ["H158", "F184"]

    elif max_z == 1.7:
        if spec_SNR == None:
            filters = ["Y106", "J129", "H158", "F184"]
        else:
            filters = ["Y106", "J129"]
            parfilters = ["H158", "F184"]

    return filters, parfilters



def make_param_file(flname, spec_SNR, square_degrees, max_z, spectra_depth_and_phase, h1rg, survey_time, TTel, include_f184_for_08):
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
    f.write("#SBATCH --mem=%i00\n" % (35))
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
    f.write("total_survey_time,%.2f\n" % survey_time)
    f.write("survey_duration,2.5\n") # Separate from how long the survey actually takes. This is the maximum!
    f.write("hours_per_visit,15.\n")
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


    f.write("__________________________________________\n")

    f.write("tier_name,Tier\n")
    f.write("tier_fraction_time,1.0\n")
    f.write("square_degrees,%.2f\n" % square_degrees)


    filters, parfilters = get_filters(include_f184_for_08, max_z, spec_SNR)

    f.write("parallel_filters,%s\n" % (",".join(parfilters)))




    dithers_per_filter = [1]*len(filters)
    dithers_per_filter = [str(item) for item in dithers_per_filter]
    exp_times_per_dither = [str(get_exptime(item, max_z)) for item in filters]


    f.write("filters," + ",".join(filters) + "\n")
    f.write("exp_times_per_dither," + ",".join(exp_times_per_dither) + "\n")
    f.write("cadence,5\n")
    f.write("dithers_per_filter," + ",".join(dithers_per_filter) + "\n")

    if spec_SNR == None:
        trigger_redshifts = [0, 0.01, 3]
        trigger_fraction = [1, 0, 0]
    else:
        if max_z < 0.81:
            trigger_redshifts = [0, max_z, max_z+0.01, 3]
            trigger_fraction = [1, 1, 0, 0]
        else:
            trigger_redshifts = [0, 0.8, 0.81, max_z, max_z+0.01, 3]
            trigger_fraction = [1, 1, 0.5, 0.5, 0, 0]

    trigger_redshifts = ["%.2f" % item for item in trigger_redshifts]
    trigger_fraction = ["%.2f" % item for item in trigger_fraction]

    f.write("trigger_redshifts," + ",".join(trigger_redshifts) + "\n")
    f.write("trigger_fraction," + ",".join(trigger_fraction) + "\n")

        
    f.close()
    #print commands.getoutput("sbatch " + working_dir + "/survey.sh")


def lnspace(start, stop, samples):
    return list(around(exp(linspace(log(start), log(stop), samples))))

def sqspace(start, stop, samples, decimals = 0):
    return list(around(linspace(sqrt(start), sqrt(stop), samples)**2., decimals = decimals))

def get_combinations():
    params_list = []

    for max_z in [0.8, 1.2, 1.7]:
        for spec_SNR in ([3.5*sqrt(15.), 6*sqrt(15.), 10.*sqrt(15), sqrt(136.*15.)], None):
            for include_f184_for_08 in [1] + [0]*(max_z == 0.8):
                for deg_scale in [1.0] + [0.7, 0.8, 0.9]*(spec_SNR != None):
                    try:
                        spec_SNR[2]
                        spec_SNR = ["%.2f" % item for item in spec_SNR]
                    except:
                        pass
                    
                    params = dict(max_z = max_z, spec_SNR = spec_SNR, include_f184_for_08 = include_f184_for_08)
                    
                    params["spectra_depth_and_phase"] = ("shallow  -10", "medium -1", "deep 1", "shallow 5", "shallow 10")
                    
                    
                    
                    params["h1rg"] = 0
                    params["TTel"] = 260
                    
                    params["survey_time"] = 0.3
                    
                    filters, parfilters = get_filters(include_f184_for_08, max_z, spec_SNR)
                    exp_times = array([get_exptime(filt, max_z) for filt in filters + parfilters])
                    exp_times += 65.
                    
                    pointings = (18*3600.)/sum(exp_times)
                
                    square_degrees = pointings*0.28
                    if spec_SNR != None:
                        square_degrees *= deg_scale

                    params["square_degrees"] = around(square_degrees)
                    
                    print "square_degrees", params["square_degrees"], max_z, include_f184_for_08


                    if params_list.count(params) == 0:
                        params_list.append(params)

    return params_list


params_list = get_combinations()


commands.getoutput("rm -f slurm-*.out")
print "Make sure directory is empty!"



print "Ready to start..."

thecount = 1
for params in params_list:
    make_param_file(sys.argv[1] + "/survey_sqdeg=%.1f_z=%.2f_F184=%i_sp=%i/paramfile_%05i.csv" % (params["square_degrees"], params["max_z"], params["include_f184_for_08"], params["spec_SNR"] != None, thecount), **params)
    thecount += 1

