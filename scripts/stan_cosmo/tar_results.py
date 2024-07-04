import subprocess
import sys

dr = sys.argv[1]
tfl = dr.replace("/", "") + ".tar"

cmd = "tar -cvf " + tfl + " " + dr + "/*/param* " + dr + "/*/fisher_log.txt " + dr + "/*/*FoM*"
print(cmd)
print(subprocess.getoutput(cmd))
print(subprocess.getoutput("gzip " + tfl))
