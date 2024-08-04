import glob
from subprocess import getoutput
import sys

def run_job(dr, gray_disp, color_scatter_opt, color_scatter_nir, twins):
    f = open(dr + "/tmp.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --partition=shared
#SBATCH --time=0-08:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile
export WFIRST=/home/drubin/wfirst-sim/
export WFIRST_SIM_DATA=/home/drubin/wfirst-sim-data/
pip install sncosmo
pip install sep

cd """ + pwd + "/" + dr + """
python $WFIRST/scripts/stan_cosmo/STEP2_Analytic_Fisher.py survey.pickle --SNRMax 0 --model_res 9 --gray_disp """ + str(gray_disp) + " --color_scatter_opt " + str(color_scatter_opt) + " --color_scatter_nir " + str(color_scatter_nir) + " --twins " + str(twins) +  """ > fisher_log.txt
python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat*fits""")
    f.close()                                                                                                                                                                       
    print(getoutput("cd " + dr + "\nsbatch tmp.sh"))
    

pwd = getoutput("pwd")

redo = int(sys.argv[1])

comb_mats = glob.glob("*/comb_mat*fits")
# yr=0.500_wi=057_mi=023_di=020_wp=000_mp=000_dp=000_nnearby=00800_ZYJH+RZYJ+RZYJH_cad=10+05+05_PN=0_0/comb_mat__SNRMax=0_res=09_bins=025_disp=0.100_scatopt=0.030_scatnir=0.030.fits

for dr in glob.glob("*nnearby*"):
    for gray_disp, color_scatter_opt, color_scatter_nir, twins in [(0.1, 0.03, 0.03, 0), (0.08, 0.04, 0.04, 0), (0.08, 0.04, 0.02, 0), (0, 0, 0, 1)]:
        if redo:
            run_job(dr = dr, gray_disp = gray_disp, color_scatter_opt = color_scatter_opt, color_scatter_nir = color_scatter_nir, twins = twins)
        else:
            found = 0
            for comb_mat in comb_mats:
                if found == 0 and comb_mat.count(dr) == 1:
                    if twins == 0:
                        if comb_mat.count("disp=%.3f_scatopt=%.3f_scatnir=%.3f" % (gray_disp, color_scatter_opt, color_scatter_nir)) == 1:
                            found = 1
                    else:
                        if comb_mat.count("disp=1.000"):
                            found = 1
            if found == 1:
                print("Found ", (gray_disp, color_scatter_opt, color_scatter_nir, twins))
            else:
                run_job(dr = dr, gray_disp = gray_disp, color_scatter_opt = color_scatter_opt, color_scatter_nir = color_scatter_nir, twins = twins)

