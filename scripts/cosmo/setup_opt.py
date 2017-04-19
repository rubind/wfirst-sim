import commands

f = open("opt.sh", 'w')

binw = "'" + '"Binw"' + "'"
detf = "'" + '"DETF"' +"'"

for slew_time in [0, 70]:
    for start in [0, 1]:
        for FoM_type in ["[" + detf + "]", "[" + binw + ",0,0.15,0.5]", "[" + binw + ",0,0.1,0.2,0.75]"]:
            # One day, systematcs and typing high-redshift bins together and exploring max z

            wd = "opt_slew=%.1f_start=%i_FoM=%s" % (slew_time, start, FoM_type.replace("[", "").replace("]", "").replace(",", "+").replace("'", "").replace('"', ""))

            commands.getoutput("rm -fr " + wd)
            commands.getoutput("mkdir " + wd)
            commands.getoutput("cp paramfile_wrap.txt " + wd)

            f.write("cd " + commands.getoutput("pwd") + "/" + wd + '\n')
            f.write("~/anaconda2/bin/python ../optimize_survey.py %.2f %i %s > log.txt & \n" % (slew_time, start, FoM_type))
f.write("wait")

