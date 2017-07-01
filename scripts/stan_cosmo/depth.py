from numpy import *
import cPickle as pickle
import matplotlib.pyplot as plt
import commands
from scipy.interpolate import interp1d
import sys

SN_data = pickle.load(open(sys.argv[1], 'rb'))


for test_tier in SN_data["test_points"]:
    for key in test_tier:
        print key
    ra_dec_set = zip(test_tier["RAs"], test_tier["Decs"])
    ra_dec_set = list(set(ra_dec_set))
    print ra_dec_set

    print len(test_tier["Observations"])
    
    SNR_by_pos = [{} for i in range(len(ra_dec_set))]
    mag_by_pos = [{} for i in range(len(ra_dec_set))]
    
    for i in range(len(test_tier["Observations"])):
        ind = ra_dec_set.index((test_tier["RAs"][i], test_tier["Decs"][i]))
        
        for filt in test_tier["Observations"][i]:
            this_SNR = sqrt(float(
                            dot(test_tier["Observations"][i][filt], test_tier["Observations"][i][filt])
                            ))
            this_mag = test_tier["Mags"][i]
                
            if SNR_by_pos[ind].has_key(filt):
                SNR_by_pos[ind][filt] = append(SNR_by_pos[ind][filt], this_SNR)
                mag_by_pos[ind][filt] = append(mag_by_pos[ind][filt], this_mag)
            else:
                SNR_by_pos[ind][filt] = array([this_SNR])
                mag_by_pos[ind][filt] = array([this_mag])

    depths_by_filt = {}

    for i in range(len(SNR_by_pos)):
        for filt in SNR_by_pos[i]:
            if any(SNR_by_pos[i][filt] > 0):
                inds = argsort(SNR_by_pos[i][filt])
                ifn = interp1d(SNR_by_pos[i][filt][inds], mag_by_pos[i][filt][inds], kind = 'linear')
                
                try:
                    depth = ifn(5.)
                except:
                    print "Couldn't get depth!", SNR_by_pos[i][filt][inds]
                    depth = max(mag_by_pos[i][filt])
                    
                if depths_by_filt.has_key(filt):
                    depths_by_filt[filt].append(depth)
                else:
                    depths_by_filt[filt] = [depth]
    for filt in depths_by_filt:
        print filt, "%.1f" % median(depths_by_filt[filt])
