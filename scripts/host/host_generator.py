from numpy import *
from scipy.interpolate import interp1d
import os
from scipy.ndimage.filters import gaussian_filter1d

wfirst_data_path = os.environ["WFIRST_SIM_DATA"]

def get_ST_flux(filter_fn, templ_fn, redshift):
    waves = arange(3000., 18000., 5.)
    
    flux = sum(filter_fn(waves)*waves*templ_fn(waves/(1. + redshift)))
    stflux = sum(filter_fn(waves)*waves)
    
    return flux/stflux


def make_galaxy_spectrum(redshifts, median_not_random = False, show_plots = False, convolve = True):
    data = loadtxt(wfirst_data_path + "/host/templ.txt")
    wave = data[:,0]
    templ = data[:,1]
    onecomp = data[:,2]


    f160w_ST_mag_per_arcsec = 24.69330219319535 + random.normal(size = len(redshifts))*2.1544518737582403*(1 - median_not_random)
    f160w_flamb_per_arcsec = 10.**(-0.4*(21.1 + f160w_ST_mag_per_arcsec))

    obswaves = exp(arange(log(3900.), log(25500.), 0.0005))
    

    templs = []
    for i in range(len(redshifts)):
        temptempl = templ*exp(random.normal()*onecomp*(1 - median_not_random))
        temptempl = interp1d(wave*(1 + redshifts[i]), temptempl, kind = 'linear')(obswaves)
        if convolve:
            temptempl = gaussian_filter1d(temptempl, 4.)

        ind = argmin(abs(obswaves - 15500.))
        temptempl *= f160w_flamb_per_arcsec[i]/temptempl[ind]

        if show_plots:
            inds = where((obswaves >= 4200.)*(obswaves <= 25000.))
            plt.plot(obswaves[inds], temptempl[inds], color = ['b','g','c','r','m','y'][i % 6])
        templs.append(interp1d(obswaves[::4], temptempl[::4], kind = 'linear'))
            
    if show_plots:
        zodi = loadtxt(wfirst_data_path + "/pixel-level/input/aldering.txt")
        zodifn = interp1d(zodi[:,0], zodi[:,1], kind = 'linear')
        plt.plot(obswaves[inds], 10.**(zodifn(obswaves[inds])), color = 'k', linewidth = 2)
        plt.text(10000., 10.**(zodifn(10000.))*1.4, "Zodiacal", color = 'k', rotation = -8.5, bbox=dict(facecolor='w', edgecolor = 'w'))
        

        plt.xlabel("Wavelength ($\\AA$)")
        plt.ylabel("$f_{\lambda}$ per Square Arcsecond")
        plt.yscale('log')
        plt.savefig("templs.eps", bbox_inches = 'tight')
        plt.close()
    return templs


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "serif"

    make_galaxy_spectrum(random.random(size = 20)*2., show_plots = True)
