from subprocess import getoutput
import numpy as np

grepout = getoutput("grep FoM_0.26 */FoM*disp=0.100*").split('\n')
FoM = [float(item.split(None)[-1]) for item in grepout]


inds = np.argsort(FoM)



for ind in inds:
    parsed = grepout[ind].replace("_", " ").replace("=", " ").split(None)
    #print(parsed)
    

    #['yr', '0.500', 'w', '065', 'm', '005', 'd', '030', 'wp', '016', 'mp', '000', 'dp', '035', 'nnearby', '00800', 'RZJ+RZYJHF+ZYJH', 'cad', '10+05+05', 'PN', '0', '0/FoM', 'model', 'res', '09.txt:FoM', '0.26', 'comb', 'mat.fits', '288.7239454040043']

    print("years: " + parsed[1] + ", % wide/medium/deep: " +  parsed[3] + "/" + parsed[5] + "/" + parsed[7] + ", % with prism wide/medium/deep: " + parsed[9] + "/" + parsed[11] + "/" + parsed[13] +
          ", filters: " + parsed[16].replace("+", "/") + ", cadence: " + parsed[18].replace("+", "/") + ", FoM: %.1f" %  float(parsed[-1]))
