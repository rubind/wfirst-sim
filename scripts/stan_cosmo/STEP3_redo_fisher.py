import glob
from subprocess import getoutput

pwd = getoutput("pwd")

color_scatter_opt = 0.04
color_scatter_nir = 0.04


for dr in glob.glob("*nnearby*"):
    for gray_disp in [0.08]:#, 0.1]:
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
python $WFIRST/scripts/stan_cosmo/STEP2_Analytic_Fisher.py survey.pickle --SNRMax 0 --model_res 9 --gray_disp """ + str(gray_disp) + " --color_scatter_opt " + str(color_scatter_opt) + " --color_scatter_nir " + str(color_scatter_nir) +  """ > fisher_log.txt
python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat*fits""")
        f.close()

        print(getoutput("cd " + dr + "\nsbatch tmp.sh"))
    
