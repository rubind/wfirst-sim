import subprocess
import numpy as np

def run_job(train, zpunc, slopeunc, crnlunc, waveunc, wd):
    subprocess.getoutput("mkdir " + wd)
    
    f = open(wd + "/step.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --partition=kill-shared
#SBATCH --time=0-05:00:00 ## time format is DD-HH:MM:SS
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
cd """ + pwd + "/" + wd + '\n')

    f.write("python $WFIRST/scripts/stan_cosmo/STEP2_Analytic_Fisher.py ../survey.pickle --SNRMax 0 --model_res 9  --train 1 --calib 1 --color_scatter_opt 0.04 --color_scatter_nir 0.02 --gray_disp 0.08 --train %i --zpunc %f --slopeunc %f --crnlunc %f --waveunc %f \n"
            % (train, zpunc, slopeunc, crnlunc, waveunc))
    f.close()
    
    print(subprocess.getoutput("sbatch " + wd + "/step.sh"))


    
pwd = subprocess.getoutput("pwd")

for train in [0, 1]:
    for zpunc in [0.0001, 0.005]:
        for slopeunc in [0.0001, 0.007]:
            for crnlunc in [0.00000125, 0.000125]:
                for waveunc in [0.01, 0.5]:
                    run_job(train = train, zpunc = zpunc, slopeunc = slopeunc, crnlunc = crnlunc, waveunc = waveunc, wd = "step_by_step")
    


for train in [0, 1]:
    for zpunc in np.linspace(0.0001, 0.01, 20):
        run_job(train = train, zpunc = zpunc, slopeunc = 0.007, crnlunc = 0.000125, waveunc = 0.5, wd = "zpunc")
    
    for slopeunc in np.linspace(0.0001, 0.02, 20):
        run_job(train = train, zpunc = 0.005, slopeunc = slopeunc, crnlunc = 0.000125, waveunc = 0.5, wd = "slopeunc")

    for	crnlunc in np.linspace(0.00000125, 0.0005, 20):
        run_job(train = train, zpunc = 0.005, slopeunc = 0.007, crnlunc = crnlunc, waveunc = 0.5, wd = "crnlunc")

    for waveunc in np.linspace(0.01, 2., 20):
        run_job(train = train, zpunc = 0.005, slopeunc = 0.007, crnlunc = 0.000125, waveunc = waveunc, wd = "waveunc")

