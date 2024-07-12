import glob
try:
    from commands import getoutput
except:
    from subprocess import getoutput

pwd = getoutput("pwd")
    
for pfl in glob.glob("*/*pickle"):
    dr = pfl.split("/")[0]
    if glob.glob(dr + "/FoM.txt") == []:
        print("Redo ", dr)
        

        f = open(dr + "/redo.sh", 'w')
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
python $WFIRST/scripts/stan_cosmo/STEP2_Analytic_Fisher.py survey.pickle --SNRMax 0 > fisher_log.txt
python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat.fits > FoM.txt
python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat_no_model.fits > FoM_no_model.txt
python $WFIRST/scripts/stan_cosmo/FoM.py comb_mat_stat_only.fits > FoM_stat_only.txt
""")
        f.close()

        getoutput("cd " + dr + "\nsbatch redo.sh")
