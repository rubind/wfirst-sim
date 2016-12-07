import commands
import sys
from numpy import *
import os
wfirst_path = os.environ["WFIRST"]


def permute(**kwargs):
    key_list = kwargs.keys()
    print key_list
    assert len(key_list) > 1

    
    return_items = [{}]
    
    for key in key_list:
        return_items[-1][key] = kwargs[key][0]

    print "default configuration:", return_items[-1]

    for key1 in key_list:
        for item in kwargs[key1][1:]:
            return_items.append({})
            return_items[-1][key1] = item

            for key2 in key_list:
                if key1 != key2:
                    return_items[-1][key2] = kwargs[key2][0]
    print len(return_items), "permutations!"

    return_names = []
    for item in return_items:
        return_names.append("_".join([key + "=" + str(item[key]) for key in key_list]))
    print 

    return return_items, return_names


dosub = 1
one_permute = 0


if sys.argv.count("nosub"):
    dosub = 0
    del sys.argv[sys.argv.index("nosub")]

surveys = sys.argv[1:]

if os.path.isfile(surveys[0]):
    print "Reading in file!"
    surveys = [line.rstrip('\n') for line in open(sys.argv[1])]

    if surveys[-1] == "":
        del surveys[-1]

# Start with the default

if one_permute:

    permutations, names = permute(graydisp = [0.08],   
                                  nredcoeff = [2],             
                                  IPC = [0.02],
                                  crnl = [0.003],
                                  fundcal = [0.005],
                                  TTel = [260], #260
                                  redwave = [8600.], # 12800
                                  MWZP = [0.003],
                                  MWnorm = [0.05],
                                  nnearby = [800],
                                  include_sys = [1])
else:
    permutations, names = permute(graydisp = [0.08, 0.1, 0.06],   
                                  nredcoeff = [2, 3, 4],             
                                  IPC = [0.02],
                                  crnl = [0.003, 0.001, 0.002, 0.005, 0.01],

                                  fundcal = [0.005, 0.003, 0.001, 0.002, 0.01],
                                  TTel = [260, 284], #260

                                  redwave = [8600.], # 12800
                                  MWZP = [0.003],
                                  nnearby = [400, 800, 1200],
                                  MWnorm = [0.05],                                                            
                                  include_sys = [1, 0])


for permutation, name in zip(permutations, names):
    print name
    print permutation


for permutation, name in zip(permutations, names):
    # Start with default configuration

    for survey in surveys:

        wd = survey + "/FoM_" + name
        print wd

        commands.getoutput("mkdir -p " + wd)

        f = open(wd + "/unity.sh", 'w')
        f.write("#!/bin/bash -l\n")

        f.write("#SBATCH --partition=shared\n")
        f.write("#SBATCH -n 8\n")
        #f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --time=48:00:00\n")
        if commands.getoutput("hostname").count("cori"):
            f.write("#SBATCH -C haswell\n")
        f.write("#SBATCH --job-name=unity\n")
        f.write("#SBATCH --mem=15000\n")

        sys_scale = clip(permutation["include_sys"], 0.01, 1)
        nrestlamb = int(around(log(permutation["redwave"]/3300.)*21.))
        print "nrestlamb ", nrestlamb
        

        f.write("module load python/2.7-anaconda\n")
        f.write("cd " + commands.getoutput("pwd") + "/" + wd + "\n")
        f.write("srun -n 1 -c 8 python " + wfirst_path + "/scripts/stan_cosmo/STEP2_UNITY.py -p ../pickle*t -nrestlamb " + str(nrestlamb) + " -neigen 1 "
                + " -gray " + str(permutation["graydisp"])
                + " -nredcoeff " + str(permutation["nredcoeff"])
                + " -IFCIPC " + str(permutation["IPC"])
                + " -crnl " + str(permutation["crnl"]*sys_scale)
                + " -fund " + str(permutation["fundcal"]*sys_scale)
                + " -TTel " + str(permutation["TTel"]) + " -IFCmaxwave 21000 "*(permutation["TTel"] < 280)
                + " -mwnorm " + str(permutation["MWnorm"]*sys_scale) + " -mwZP " + str(permutation["MWZP"]*sys_scale) + " -mwRV " + str(0.2*sys_scale) + " -IGext " + str(0.25*sys_scale) + " -redwave " + str(permutation["redwave"])
                + " > log2.txt\n")
        # + " -IFCdark " + str(permutation["dark_current"])
        # + " -IFCRNfloor " + str(permutation["read_noise_floor"])


        f.close()

        if dosub:
            print commands.getoutput("sbatch " + wd + "/unity.sh")

