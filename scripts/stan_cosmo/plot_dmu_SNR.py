from FileRead import read_param
import glob
import matplotlib.pyplot as plt
import numpy as np

tiers = ["Nearby", "RubinDDF", "Wide", "Deep"]

for i, tier in enumerate(tiers):
    all_dmu = []
    all_SNR = []
    
    for resfl in glob.glob(tier + "/*/result_salt2.dat"):
        dmu = read_param(resfl, "dmu_estimate")
        if dmu != None and dmu < 1:
            lightfile = resfl.replace("result_salt2.dat", "lightfile")
            TOTALSNR = read_param(lightfile, "TOTALSNR")

            all_dmu.append(dmu)
            all_SNR.append(TOTALSNR)
    plt.plot(all_SNR, all_dmu, '.^v*'[i], color = 'mbgr'[i], label = tier)
plt.legend(loc = 'best')
plt.savefig("dmu_vs_SNR.pdf", bbox_inches = 'tight')
plt.close()

