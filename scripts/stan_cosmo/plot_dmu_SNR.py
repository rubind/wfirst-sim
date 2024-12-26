from FileRead import read_param
import glob
import matplotlib.pyplot as plt
import numpy as np
import tqdm

tiers = ["Nearby", "RubinDDF", "Deep", "Wide"]

for i, tier in enumerate(tiers):
    all_dmu = []
    all_SNR = []
    
    for resfl in tqdm.tqdm(glob.glob(tier + "/*/result_salt2.dat")):
        dmu = read_param(resfl, "dmu_estimate")
        if dmu != None and dmu < 1:
            lightfile = resfl.replace("result_salt2.dat", "lightfile")
            TOTALSNR = read_param(lightfile, "TOTALSNR")

            all_dmu.append(dmu)
            all_SNR.append(TOTALSNR)
    plt.plot(all_SNR, all_dmu, '.^v*'[i], color = 'mbgr'[i], label = tier, zorder = -len(all_SNR))

plt.axvline(40, color = 'k')
plt.text(50, 0.7, "Passes Minimum S/N $\\rightarrow$", )
plt.xlabel("Quadrature Sum of Light-Curve S/N")
plt.ylabel("Distance Uncertainty from SALT3 NIR")
plt.xscale('log')
plt.yscale('log')
plt.xlim(10, plt.xlim()[1])
plt.ylim(plt.ylim()[0], 1)

plt.legend(loc = 'best')
plt.savefig("dmu_vs_SNR.pdf", bbox_inches = 'tight')
plt.close()

