from scipy.interpolate import RectBivariateSpline, interp1d
from numpy import *
from DavidsNM import save_img
from astropy import cosmology
from FileRead import readcol
from extinction import fm07, ccm89
import os

def NMAD(vals, axis = None):
    return 1.4826*median(abs(vals - median(vals, axis = axis)), axis = axis)

def get_weights(true_x, projs_0, dP, NMAD_projs, clip_kernel):
    weights = exp(-0.5* sum((true_x - projs_0)**2. /(dP*NMAD_projs)**2., axis = 1))
    if clip_kernel:
        weights *= (weights > 0.0001)
    return weights


def file_to_fn(fl):
    [x, y] = readcol(fl, 'ff')
    return interp1d(x, y, kind = 'linear', fill_value = 0, bounds_error = False)

def get_data(dP = 0.5, clip_kernel = 0):

    wfirst_path = os.environ["WFIRST"]

    [phase, wave] = readcol(wfirst_path + "/../ClareModel/sample_sn.dat", 'ff')
    
    data = loadtxt(wfirst_path + "/../ClareModel/EMfa_eigenvectors_Aprem_15.dat")
    projs = loadtxt(wfirst_path + "/../ClareModel/EMfa_projections_Aprem_15.dat")

    projs[:, 0] *= -1
    data[:, 0] *= -1



    ncomp = len(projs[0])

    for i in range(1, ncomp):
        projs[:, i] = projs[:, i]/projs[:, 0]
    

    phases = sort(unique(phase))
    waves = sort(unique(wave))
    

    print "data.shape, len(phases), len(waves)", data.shape, len(phases), len(waves)
    print "projs.shape", projs.shape

    compfns = []
    for i in range(len(data[0])):
        datares = reshape(data[:,i], [len(waves), len(phases)])
        save_img(datares, str(i) + ".fits")
        compfns.append(RectBivariateSpline(waves, phases, datares))

    print phases

    norms = projs[:, 0]
    projs = projs[:, 1:]

    NMAD_projs = NMAD(projs, axis = 0)

    norms_in_mags = -2.5*log10(norms)

    smoothed_norms_in_mags = []
    for i in range(len(norms)):
        weights = get_weights(projs[i], projs, dP = dP, NMAD_projs = NMAD_projs, clip_kernel = clip_kernel)
        smoothed_norms_in_mags.append(sum(norms_in_mags*weights)/sum(weights))
    smoothed_norms_in_mags = array(smoothed_norms_in_mags)
    smoothed_norms = 10.**(-0.4*smoothed_norms_in_mags)

    return compfns[0], norms, compfns[1:], projs, smoothed_norms, NMAD_projs



def get_filts():
    filtfns = {}
    filtwaves = {}
    test_waves = arange(3000., 25000.)

    wfirst_data_path = os.environ["WFIRST_SIM_DATA"]


    for filt in ["g", "r", "i", "z", "Y",
                 "R062", "Z087", "Y106", "J129", "H158", "F184"]:
        if len(filt) > 1:
            filtfns[filt] = file_to_fn(wfirst_data_path + "/pixel-level/input/" + filt + ".txt")
        else:
            filtfns[filt] = file_to_fn(wfirst_data_path + "/pixel-level/input/LSST_" + filt + ".txt")

        ifilt = filtfns[filt](test_waves)
        ifilt = ifilt[1:]*sign(ifilt[1:] - ifilt[:-1])
        print filt
        ind = argmin(abs(ifilt - ifilt.max()/3.))
        minwave = test_waves[1:][ind]
        ind = argmin(abs(ifilt + ifilt.max()/3.))
        maxwave = test_waves[1:][ind]
        filtwaves[filt] = array([minwave, maxwave])

    filtnames = {}
    for filt in filtwaves:
        if filt.count("LSST"):
            filtnames[filt] = "LSST::" + filt.split("_")[-1]
        else:
            filtnames[filt] = "WFIRST::" + filt

    return filtfns, filtwaves, filtnames


def get_AB_mags(true_norm, true_proj, true_AV, z, filtfn, phases, lambs, ncomp, meanfn, compfns, cosmo = None):
    if cosmo == None:
        cosmo = cosmology.FlatLambdaCDM(Om0 = 0.3, H0 = 70.)

    mu = cosmo.distmod(z).value
    #print "mu ", z, mu, 10.**(-0.4*mu)
    
    tmp_restlambs = lambs/(1. + z)
    tmp_nphase = len(phases)

    assert len(true_proj) == ncomp
    flux = true_norm*meanfn(tmp_restlambs, [phases])
    for j in range(ncomp):
        flux += true_norm*true_proj[j]*compfns[j](tmp_restlambs, [phases])
    flux *= outer(10.**(-0.4*ccm89(a_v = true_AV, r_v = 3.1, wave = tmp_restlambs)) / tmp_restlambs, ones(tmp_nphase)) # Divide by tmp_restlambs for normalization
    flux /= (1. + z)
    
    
    flux *= 10.**(-0.4*mu)
    

    #return array([-2.5*log10(
    return array([(flux[:,i]*lambs*filtfn(lambs)).sum()/(filtfn(lambs)/lambs).sum()
                  for i in range(len(phases))])
    #) for i in range(len(phases))])


def get_spectrum(true_norm, true_proj, true_AV, restlambs, ncomp, meanfn, compfns):
    flux = true_norm*meanfn(restlambs, 0)

    for j in range(ncomp):
        flux += true_norm*true_proj[j]*compfns[j](restlambs, 0)
    flux = flux[:,0]

    flux *= 10.**(-0.4*ccm89(a_v = true_AV, r_v = 3.1, wave = restlambs)) 
    return flux

