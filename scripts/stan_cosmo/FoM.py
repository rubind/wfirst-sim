import numpy as np
from astropy.io import fits
import os
import sys
import matplotlib.pyplot as plt
import tqdm


wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/cosmo/")
from astro_functions import get_FoM


fgrid = open("FoM_grid.txt", 'w')

for item in tqdm.tqdm(sys.argv[1:]):
    try:
        f = fits.open(item)
        combmat = f[0].data
        print(combmat.shape)
        cmat = combmat[1:]
        z_list = combmat[0]
        
        f.close()
        bad_file = 0
    except:
        bad_file = 1
        print("bad_file", item)

    if bad_file == 0:
        print("z_list", z_list)
        
        
        
        
        f = open(item.replace(".fits", ".txt").replace("comb_mat", "FoM_comb_mat"), 'w')
        FoM = get_FoM(cmat, z_list, shift_constraint = 0.0035)[0]
        f.write("FoM_0.35 " + str(item) + " " + str(FoM) + '\n')
        
        FoM = get_FoM(cmat, z_list, shift_constraint = 0.0026)[0]
        print("FoM_0.26 ", item, FoM)
        f.write("FoM_0.26 " + str(item) + " " + str(FoM) + '\n')
        fgrid.write("FoM_0.26 " + str(item) + " " + str(FoM) + '\n')
        
        if len(sys.argv) < 20:
            plt.plot(z_list, np.sqrt(np.diag(cmat)), label = item + "\nFoM_0.26: %.1f" % FoM)

    
        FoM = get_FoM(cmat, z_list, shift_constraint = 0.002)[0]
        print("FoM_0.2 ", item, FoM)
        f.write("FoM_0.2 " + str(item) + " " + str(FoM) + '\n')
        f.close()

fgrid.close()
    
if len(sys.argv) < 20:
    plt.legend(loc = 'best')
    plt.savefig("dmu_vs_z.pdf", bbox_inches = 'tight')
    plt.close()
