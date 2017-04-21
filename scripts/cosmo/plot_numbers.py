import commands
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import glob
from numpy import *


drs = glob.glob("opt_slew*run=1_*")
drs.sort()


plt.figure(figsize = (24, 10))
for i in range(len(drs)):
    plt.subplot(2, len(drs)/2,
                i+1)

    for start in [1, 2]:

        for pos in range(-4, 0, 1):
            print drs[i]
            grepout = commands.getoutput("grep scaled " + drs[i].replace("run=1", "run=%i" % start) + "/log.txt | grep '\[' ").split('\n')[pos]

            
            grepout = grepout.replace("scaled", "")
        
            try:
                nsne = eval(grepout)
                print nsne
            except:
                pass

        for j in range(len(nsne)):
            plt.fill_between([j*0.1 + 0.1, (j + 1)*0.1 + 0.1], [0,0], [nsne[j]]*2, color = 'grb'[start-1], alpha = 0.5, label = (j == 0)*"")
        plt.title(drs[i])

plt.savefig("NSNe.pdf")
plt.close()
