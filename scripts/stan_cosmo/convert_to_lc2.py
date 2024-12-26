import pickle
import numpy as np
import subprocess
import tqdm
import sys


def write_header(f_lcs, filt):
    if len(filt) == 1:
        instrument = "LSST"
        band = "LSST_" + filt
    else:
        instrument = "WFIRST"
        band = filt
    
    f_lcs[filt] = open(wd + "lc2fit_" + filt + ".dat", 'w')
    f_lcs[filt].write("""#Date :
#Flux :
#Fluxerr :
#ZP :
#end :
@INSTRUMENT  """ + instrument + """
@BAND  """ + band + """
@MAGSYS  AB
""")
    return f_lcs


frac_to_fit = float(sys.argv[1])

SN_data = pickle.load(open(sys.argv[2], 'rb'))

for key in SN_data:
    print("SN_data", key)
    
print("SN_observations[0]", SN_data["SN_observations"][0])
print("SN_table", SN_data["SN_table"])


assert len(SN_data["SN_table"]["survey_fields"]) == len(SN_data["SN_observations"])

subprocess.getoutput("rm -fr LCs")



for i in tqdm.trange(len(SN_data["SN_observations"])):
    if np.random.random() < frac_to_fit:
        this_obs = SN_data["SN_observations"][i]

        SNRs = this_obs["fluxes"]/this_obs["dfluxes"]
        inds = np.where(SNRs > 2)
        SNRs = SNRs[inds]
        TOTALSNR = np.sqrt(np.dot(SNRs, SNRs))

        if len(SNRs) > 3 and sum(SNRs**2.) > 10**2.:
            wd = "LCs/" + SN_data["SN_table"]["survey_fields"][i] + "/SN_%05i/" % i
            subprocess.getoutput("mkdir -p " + wd)


            f = open(wd + "lightfile", 'w')
            f.write("MWEBV 0.0\n")
            f.write("Mass " + str(this_obs["host_mass"]) + "  -0.1  0.1\n")
            f.write("RA  0.0 \n")
            f.write("Dec  0.0 \n")
            f.write("z_heliocentric  " + str(SN_data["SN_table"]["redshifts"][i]) + '\n')
            f.write("z_CMB  " + str(SN_data["SN_table"]["redshifts"][i]) + '\n')
            f.write("TOTALSNR  " + str(TOTALSNR) + '\n')
            f.close()
            
            f_lcs = {}

        
            for j in range(len(this_obs["dates"])):
                filt = this_obs["filts"][j]
                if filt not in f_lcs:
                    f_lcs = write_header(f_lcs, filt)

                f_lcs[filt].write("  ".join([str(this_obs["dates"][j]),
                                             str(this_obs["fluxes"][j]),
                                             str(this_obs["dfluxes"][j]),
                                             "25.0"]) + '\n')
            for filt in f_lcs:
                f_lcs[filt].close()
        
        
f = open("LCs/fit.sh", 'w')
f.write("python $PATHMODEL/python_code/slurmfit.py 10 dontsort -v salt3nir_final -w 3000 16000")
f.close()
