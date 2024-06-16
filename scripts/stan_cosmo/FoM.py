import numpy as np
from astropy.io import fits
import os
import sys
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/cosmo/")
from astro_functions import get_FoM

for item in sys.argv[1:]:
    f = fits.open(item)
    combmat = f[0].data
    print(combmat.shape)
    cmat = combmat[1:]
    z_list = combmat[0]
    
    f.close()

    print("z_list", z_list)



    FoM = get_FoM(cmat, z_list, shift_constraint = 0.0026)[0]
    print("FoM_0.26 ", item, FoM)
    
    FoM = get_FoM(cmat, z_list, shift_constraint = 0.002)[0]
    print("FoM_0.2 ", item, FoM)
        
