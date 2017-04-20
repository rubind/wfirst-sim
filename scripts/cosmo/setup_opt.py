import commands



binw = "'" + '"Binw"' + "'"
detf = "'" + '"DETF"' +"'"

for slew_time in [0, 70]:
    for run in range(3):
        for FoM_type in ["[" + detf + "]", "[" + binw + ",0,0.15,0.5]", "[" + binw + ",0,0.1,0.2,0.75]"]:
            # One day, systematcs and typing high-redshift bins together and exploring max z

            wd = "opt_slew=%.1f_run=%i_FoM=%s" % (slew_time, run+1, FoM_type.replace("[", "").replace("]", "").replace(",", "+").replace("'", "").replace('"', ""))

            commands.getoutput("rm -fr " + wd)
            commands.getoutput("mkdir " + wd)
            commands.getoutput("cp paramfile_wrap.txt " + wd)

            f = open("tmp.sh", 'w')
            f.write("""#!/bin/bash -l
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --job-name=opt
#SBATCH --mem=2500
module load python/2.7-anaconda
cd """ + commands.getoutput("pwd") + "/" + wd + """
srun -n 1 -c 1 python ../optimize_survey.py """ + "%.2f %s > log.txt & \n" % (slew_time, FoM_type))
            f.close()

            print commands.getoutput("sbatch tmp.sh")

