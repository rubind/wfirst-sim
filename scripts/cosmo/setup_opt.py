import commands
import sys

def liststr(the_list):
    return str(the_list).replace(",", "").replace("[", "").replace("]", "")

hostname = commands.getoutput("hostname")


binw = "'" + '"Binw"' + "'"
binrho = "'" + '"Binrho"' + "'"
detf = "'" + '"DETF"' +"'"

together = int(sys.argv[1])
parallel = int(sys.argv[2])

if together:
    f = open("tmp.sh", 'w')

runs_written = 1

for slew_time in [70.]:
    for LSST in [0, 5]:
        for run in range(2):
            for FoM_type in [("DETF", []), ("Binw", [0,0.15,0.5]), ("Binw", [0,0.1,0.2,0.75]), ("Binrho", [0,0.25,0.5,1.0,2.0])]:
                for maxz in [1.8] + [1.5, 1.2, 0.9]*(FoM_type[0] == "DETF"):
                # One day, systematcs
                
                    wd = "opt_slew=%.1f_run=%i_FoM=%s_LSST=%i_maxz=%.2f" % (slew_time, run+1, str(FoM_type).replace("[", "").replace("]", "").replace(",", "+").replace("'", "").replace('"', "").replace(" ", "").replace("(", "").replace(")", ""), LSST, maxz)

                    commands.getoutput("rm -fr " + wd)
                    commands.getoutput("mkdir " + wd)
                    commands.getoutput("cp paramfile_wrap.txt " + wd)

                    if not together:
                        f = open("tmp.sh", 'w')

                    if hostname.count("rivoli"):
                        f.write("cd " + commands.getoutput("pwd") + "/" + wd + '\n')
                        f.write("unbuffer ~/anaconda2/bin/python ../optimize_survey.py " +
                                "-slew %.2f -FoM %s -bins %s -LSST %i -maxz %.2f > log.txt %s \n" % (slew_time, FoM_type[0], liststr(FoM_type[1]), LSST, maxz, "&"*parallel))

                        if not together:
                            f.close()
                            print commands.getoutput("qsub tmp.sh")

                    else:
                        if together:
                            f.write("cd " + commands.getoutput("pwd") + "/" + wd + '\n')
                            f.write("unbuffer " + "~/anaconda2/bin/"*(1 - hostname.count("Davids")) +
                                    "python ../optimize_survey.py " + "-slew %.2f -FoM %s -bins %s -LSST %i -maxz %.2f  > log.txt %s \n" % (slew_time, FoM_type[0], liststr(FoM_type[1]), LSST, maxz, "&"*parallel))
                            f.write("cd ..\n")
                            if runs_written % 2 == 0:
                                #f.write("wait\n")
                                pass
                            runs_written += 1
                        else:
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

if together:
    f.close()
