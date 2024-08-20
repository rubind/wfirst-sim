from subprocess import getoutput
import numpy as np
import sys



#grepout = getoutput("grep FoM_0.26 */FoM*disp=0.100*").split('\n')

grepout = getoutput("grep FoM_0.26 " + sys.argv[1] + " | grep 'scatnir\=0.02'").replace("FoM_0.26", "").replace("e+", "e") # + " | grep scatter_nir\=0.02").split('\n')

grepout = grepout.split('\n')
#print(grepout[:4])

FoM = [float(item.split(None)[-1]) for item in grepout]


inds = np.argsort(FoM)


parsed = grepout[0].replace("_", " ").replace("=", " ").split(None)
for i, x in enumerate(parsed):
    print(i, x)


for ind in inds:
    parsed = grepout[ind].replace("_", " ").replace("=", " ").split(None)
    #print(parsed)
    

    #['yr', '0.500', 'w', '065', 'm', '005', 'd', '030', 'wp', '016', 'mp', '000', 'dp', '035', 'nnearby', '00800', 'RZJ+RZYJHF+ZYJH', 'cad', '10+05+05', 'PN', '0', '0/FoM', 'model', 'res', '09.txt:FoM', '0.26', 'comb', 'mat.fits', '288.7239454040043']

    #for i, x in enumerate(parsed):
    #    print(i, x)
    #ffff

    tier_percents = [int(parsed[3]), int(parsed[5]), int(parsed[7])]
    filters = []
    cadences = []
    redshifts = []
    areas = []

    for i in range(3):
        if tier_percents[i] > 1:
            filters.append(parsed[16].split("+")[i])
            cadences.append(parsed[18].split("+")[i])
            redshifts.append(parsed[24].split("+")[i])
            areas.append(parsed[26].split("+")[i])
        else:
            filters.append("-")
            cadences.append("-")
            redshifts.append("-")
            areas.append("-")
    
    
    print("years: " + parsed[1] + ", % wide/medium/deep: " +  parsed[3] + "/" + parsed[5] + "/" + parsed[7] + ", % with prism wide/medium/deep: " + parsed[9] + "/" + parsed[11] + "/" + parsed[13] + (" (%02i total)" % (float(parsed[9]) + float(parsed[11]) + float(parsed[13])))  +
          ", filters: " + "/".join(filters) + ", cadence: " + "/".join(cadences) + ", targeting redshifts:", "/".join(redshifts) + ", areas: " + "/".join(areas) + ", FoM: %.1f" %  float(parsed[-1]))
