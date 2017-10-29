from numpy import *
from astropy.io import fits as pyfits
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d
import glob
import commands
import sys
from astropy import cosmology
import sncosmo
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
from matplotlib import style
#style.use('dark_background')
dark_background = False#True
import sep
import types
import time
import os
wfirst_path = os.environ["WFIRST_SIM_DATA"]


######################################### Helper FNs #########################################
def save_img(dat, imname, waves = None):

    commands.getoutput("rm -f " + imname)
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname)
    fitsobj.close()


def interpfile(the_file, norm = False, bounds_error = True):

    if the_file.count(".fits"):
        f = pyfits.open(the_file)
        data = f[0].data
        f.close()
    else:
        data = loadtxt(the_file)
        if norm:
            data[:,1] /= data[:,1].max()
    return interp1d(data[:,0], data[:,1], kind = 'linear', fill_value = 0., bounds_error = bounds_error)


def bin_vals_fixed_bins(xvals, yvals, yerrs, xbins):
    yweights = 1./yerrs**2.

    closest_bins = array([argmin(abs(xbins - item)) for item in xvals])
    assert len(closest_bins) == len(xvals), "Problem with assigning wavelengths!"
    assert len(xvals) == len(yvals), str(len(xvals)) + " " + str(len(yvals))
    assert len(xvals) == len(yerrs)

    binnedy = zeros(len(xbins), dtype=float64)
    binnedw = zeros(len(xbins), dtype=float64)
    for i in range(len(xvals)):
        binnedy[closest_bins[i]] += yvals[i]*yweights[i]
        binnedw[closest_bins[i]] += yweights[i]
    return binnedy/binnedw, 1./sqrt(binnedw)


def convert_to_flamb_per_square_arcsec(lamb_in_micron, flux_in_MJy_per_steradian):
    newlamb = lamb_in_micron*10000. # Convert to A
    newflux = flux_in_MJy_per_steradian*1.e6  # Convert to Jy
    
    newflux /= 3630.780547 # Convert to AB
    newflux /= 42545170296.2 # Convert to square arcsec
    newflux *= (0.10884806248/newlamb**2.)
    
    return newlamb, newflux

def convert_to_MJy_per_steradian(lamb_in_A, flamb_per_square_arcsec):
    newlamb = lamb_in_A/1.e4 # Convert to microns

    newflux = flamb_per_square_arcsec / (0.10884806248/lamb_in_A**2.) # Convert to fnu
    newflux *= 42545170296.2 # Convert to steradians
    newflux *= 3630.780547 # Convert to Jy
    newflux = newflux / 1.e6  # Convert to MJy

    return newlamb, newflux


def integrate_spec(waves, dwaves, redshift, restframe_bins, f_lamb_SN, yerr_flamb):
    signal_to_noises = {"rest_frame_band_S/N": {},
                        "rest_frame_quadrature_S/N": {},
                        "rest_frame_mean_S/N": {},
                        "rest_frame_band_mag": {}}


    for i in range(len(restframe_bins) - 1):

        if (waves[0] - 0.5*dwaves[0])/(1. + redshift) <= restframe_bins[i] and (waves[-1] + 0.5*dwaves[-1])/(1. + redshift) >= restframe_bins[i+1]:
            inds = where((waves/(1. + redshift) >= restframe_bins[i])*(waves/(1. + redshift) < restframe_bins[i+1]))
            nresl = 0.5*float(len(inds[0]))
            restkey = (restframe_bins[i], restframe_bins[i+1])

            this_errs = sqrt(sum((yerr_flamb[inds]*dwaves[inds]*waves[inds])**2.))
            this_photons = sum(f_lamb_SN[inds]*dwaves[inds]*waves[inds])


            signal_to_noises["rest_frame_band_S/N"][restkey] = this_photons/this_errs
            signal_to_noises["rest_frame_quadrature_S/N"][restkey] = sqrt(sum((f_lamb_SN[inds]/yerr_flamb[inds])**2.))
            signal_to_noises["rest_frame_mean_S/N"][restkey] = this_photons/this_errs/sqrt(nresl)
            signal_to_noises["rest_frame_band_mag"][restkey] = -2.5*log10(  sum(f_lamb_SN[inds]*waves[inds]*dwaves[inds])/sum(0.10884806248 *dwaves[inds] / waves[inds])  )
    return signal_to_noises

"""
waves = arange(6950., 9100., 75.)*2.2
dwaves = ones(len(waves))*75.
redshift = 1.2
restframe_bins = [7000., 9000]
f_lamb_SN = ones(len(waves)) * 1.
yerr_flamb = ones(len(waves)) * 0.5
print waves

print integrate_spec(waves, dwaves, redshift, restframe_bins, f_lamb_SN, yerr_flamb)
stop_here
"""

def resolution_to_wavelengths(source_dir, IFURfl, min_wave, max_wave, waves = None):
    if hasattr(IFURfl, '__call__'):
        spec_R = IFURfl
    else:
        spec_R = interpfile(source_dir + "/" + IFURfl)

    if any(waves == None):
        waves = array([min_wave], dtype=float64)
        while waves[-1] <= max_wave:
            waves = append(waves, waves[-1]*exp(0.5/spec_R(waves[-1]))
                       )
        if waves[-1] > max_wave:
            waves = waves[:-1]
        dwaves = waves / (2.*spec_R(waves))
    else:
        super_waves = concatenate(([waves[0] - (waves[1] - waves[0])], waves, [waves[-1] + (waves[-1] - waves[-2])]))
        super_interm = 0.5*(super_waves[1:] + super_waves[:-1])
        dwaves = super_interm[1:] - super_interm[:-1]
    return waves, dwaves



def multi_convolve(xvals, yvals, subsample, pixel_scale):
    """pixel scale is in arcseconds per pixel. xvals/yvals subsampled by subsample."""

    sigma_in_pixels = 0.42*xvals*1e-10 / 2.36 / (pixel_scale * 4.84814e-6)
    sigma_in_pixels = sqrt(sigma_in_pixels**2. + (1.5/sqrt(12.))**2.)

    assert sigma_in_pixels.min() > 0.1 # Sanity checks for right units
    assert sigma_in_pixels.max() < 5.

    sigma_steps = linspace(sigma_in_pixels.min(), sigma_in_pixels.max(), 10.)

    conv_list = []
    for sigma in sigma_steps:
        conv_list.append(gaussian_filter1d(yvals, sigma*subsample))
    
    final_conv = yvals*0.
    for i in range(len(final_conv)):
        ind = argmin(abs(sigma_steps - sigma_in_pixels[i]))
        final_conv[i] = conv_list[ind][i]

    return final_conv

def spectrum_to_matched_resolution(specx, specy, min_wave, max_wave, IFUR, pixel_scale, subsample = 10):
    assert hasattr(IFUR, '__call__'), "IFUR must be a function"

    waves_highres = resolution_to_wavelengths(source_dir = None, IFURfl = lambda x: IFUR(x)*subsample, min_wave = min_wave, max_wave = max_wave)[0]
    ifn = interp1d(specx, specy, kind = 'linear', bounds_error = False, fill_value = 0)
    spec_highres = ifn(waves_highres)
    
    spec_conv = multi_convolve(waves_highres, spec_highres, subsample = subsample, pixel_scale = pixel_scale)
    

    waves_lowres = resolution_to_wavelengths(source_dir = None, IFURfl = lambda x: IFUR(x), min_wave = min_wave, max_wave = max_wave)[0]
    return interp1d(waves_highres, spec_conv, kind = 'linear')(waves_lowres)

"""
specx = arange(4000., 20001., 10.)
specy = specx*0
specy[::25] = 1.
specy[::50] = 0.

spectrum_to_matched_resolution(specx, specy, min_wave = 4200., max_wave = 20000., IFUR = interpfile("input/IFU_R_160720.txt"), pixel_scale = 0.075)
stop_here
"""


######################################### Image Simulation #########################################

def initialize_PSFs(pixel_scales = [15], slice_scales = [20], PSF_source = "WebbPSF", path = wfirst_path + "/pixel-level/", min_wave = -1, max_wave = 1e20):
    PSFs = {"waves": []}

    pixels = {}

    files_found = glob.glob(path + "/make_" + PSF_source + "_PSFs/*5mas*fits")
    if len(files_found) == 0:
        files_found = glob.glob(path + "/" + PSF_source + "/*5mas*fits")
    assert len(files_found) > 0, "Couldn't find " + path + "/" + PSF_source

    for fl in sort(files_found):
        the_wave = float(fl.split("/")[-1].split("_")[0])

        if the_wave > min_wave and the_wave < max_wave:
        
            f = pyfits.open(fl)
            TTPSF = f[0].data
            assert isclose(TTPSF[len(TTPSF)/2, len(TTPSF[0])/2], TTPSF.max(), rtol = 1e-2), "Max isn't in the middle! " + str(TTPSF.shape)
            f.close()

            TTPSF = fft.fft2(TTPSF)

            for slice_scale, pixel_scale in zip(slice_scales, pixel_scales):
                if not pixels.has_key((pixel_scale, slice_scale)):
                    pixels[(pixel_scale, slice_scale)] = zeros(TTPSF.shape, dtype=float64)
                    pixels[(pixel_scale, slice_scale)][:int(ceil(pixel_scale/2.)), :int(ceil(slice_scale/2.))] += 1
                    pixels[(pixel_scale, slice_scale)][-int(floor(pixel_scale/2.)):, :int(ceil(slice_scale/2.))] += 1
                    pixels[(pixel_scale, slice_scale)][:int(ceil(pixel_scale/2.)), -int(floor(slice_scale/2.)):] += 1
                    pixels[(pixel_scale, slice_scale)][-int(floor(pixel_scale/2.)):, -int(floor(slice_scale/2.)):] += 1

                    assert pixels[(pixel_scale, slice_scale)].sum() == pixel_scale*slice_scale, "Wrong indexing! " + str((pixel_scale, slice_scale)) + " " + str(pixels[(pixel_scale, slice_scale)].sum())

                    pixels[(pixel_scale, slice_scale)] = fft.fft2(pixels[(pixel_scale, slice_scale)])


            PSFs["waves"].append(the_wave
                                 )

            for slice_scale, pixel_scale in zip(slice_scales, pixel_scales):
                PSFs[(PSFs["waves"][-1], pixel_scale, slice_scale)] = array(real(fft.ifft2(TTPSF * pixels[(pixel_scale, slice_scale)])), dtype=float64)
                PSF_sep = sep.extract(PSFs[(PSFs["waves"][-1], pixel_scale, slice_scale)].copy(order='C'), thresh = PSFs[(PSFs["waves"][-1], pixel_scale, slice_scale)].max()/10.)

                ind = argmax(PSF_sep["flux"]) # Just in case diffraction gets picked up as an object!
                PSFs[(PSFs["waves"][-1], pixel_scale, slice_scale, "sigma")] = (PSF_sep["x2"][ind] * PSF_sep["y2"][ind])**0.25
                
    PSFs["waves"] = array(PSFs["waves"])
    for slice_scale, pixel_scale in zip(slice_scales, pixel_scales):
        sigmas = []
        for wave in PSFs["waves"]:
            sigmas.append(PSFs[(wave, pixel_scale, slice_scale, "sigma")])
        sigmaFN = interp1d(PSFs["waves"], sigmas, kind = 'linear')
        PSFs[(pixel_scale, slice_scale, "sigmaFN")] = sigmaFN

    PSFs["PSF_source"] = PSF_source

    return PSFs

def aint(val):
    return int(around(val))

def get_pixelized_PSF_noIPC(PSFs, pixel_scale, slice_scale, wave, offset_par, offset_perp, psfsize = 7):
    """Interpolates monochromatic, pixel-convolved PSFs to a given wavelength and takes a sample for every pixel. No IPC."""
    assert (psfsize % 2) == 1

    inds = argsort(abs(PSFs["waves"] - wave))
    nearest_waves = PSFs["waves"][inds[:2]]
    nearest_waves = sort(nearest_waves)
    nearest_weights = [(nearest_waves[1] - wave)/(nearest_waves[1] - nearest_waves[0]), (wave - nearest_waves[0])/(nearest_waves[1] - nearest_waves[0])]    
    assert abs(dot(nearest_weights, nearest_waves) - wave) < 1, "Bad interpolation!"

    len_PSF_over_two = aint(len(PSFs[nearest_waves[0], pixel_scale, slice_scale])/2)

                                                                
    PSF_at_wave = (nearest_weights[0]*PSFs[nearest_waves[0], pixel_scale, slice_scale][aint(len_PSF_over_two - floor(psfsize/2.)*pixel_scale + offset_par): aint(len_PSF_over_two + ceil(psfsize/2.)*pixel_scale + offset_par): pixel_scale,
                                                                    aint(len_PSF_over_two - floor(psfsize/2.)*slice_scale + offset_perp): aint(len_PSF_over_two + ceil(psfsize/2.)*slice_scale + offset_perp): slice_scale] + 
                   nearest_weights[1]*PSFs[nearest_waves[1], pixel_scale, slice_scale][aint(len_PSF_over_two - floor(psfsize/2.)*pixel_scale + offset_par): aint(len_PSF_over_two + ceil(psfsize/2.)*pixel_scale + offset_par): pixel_scale,
                                                                    aint(len_PSF_over_two - floor(psfsize/2.)*slice_scale + offset_perp): aint(len_PSF_over_two + ceil(psfsize/2.)*slice_scale + offset_perp): slice_scale])
    assert len(PSF_at_wave) >= psfsize, "PSF .fits is too small!"


    pixelized_PSF = PSF_at_wave
    return pixelized_PSF




    
def get_pixelized_broadband_PSF(PSFs, scale, waves, weights, IPC, offset_par = 0, offset_perp = 0):
    """This function is used for generating broadband PSFs for a given SED."""

    PSF = 0.
    total_weights = sum(weights) # Just in case weights isn't normalized

    for i in range(len(waves)):
        PSF += (weights[i]/total_weights)*get_pixelized_PSF_noIPC(PSFs, pixel_scale = scale, slice_scale = scale, wave = waves[i], offset_par = offset_par, offset_perp = offset_perp)
    

    # WFC3 IPC is nominally 0.015
    PSF_IPC = PSF*(1. - 4*IPC)
    PSF_IPC[1:] += IPC*PSF[:-1]
    PSF_IPC[:-1] += IPC*PSF[1:]

    PSF_IPC[:,1:] += IPC*PSF[:, :-1]
    PSF_IPC[:, :-1] += IPC*PSF[:, 1:]

    return PSF_IPC


def get_approximate_pixelized_broadband_PSF(PSFs, scale, waves, weights, IPC, offset_par = 0, offset_perp = 0):
    """This function is used for generating broadband PSFs for a given SED."""

    eff_lamb = sum(weights*waves)/sum(weights)

    PSF = get_pixelized_PSF_noIPC(PSFs, pixel_scale = scale, slice_scale = scale, wave = eff_lamb, offset_par = offset_par, offset_perp = offset_perp)
    
    # WFC3 IPC is nominally 0.015
    PSF_IPC = PSF*(1. - 4*IPC)
    PSF_IPC[1:] += IPC*PSF[:-1]
    PSF_IPC[:-1] += IPC*PSF[1:]

    PSF_IPC[:,1:] += IPC*PSF[:, :-1]
    PSF_IPC[:, :-1] += IPC*PSF[:, 1:]

    return PSF_IPC


def get_1D_pixelized_sliced_PSFs(PSFs, pixel_scale, slice_scale, waves, IPC, offset_par = 0, offset_perp = 0, psfsize = 7):
    """Simulates a small IFU; for a S/N calculation, that's good enough. i = parallel j = perpendicular to slice"""


    if PSFs.has_key((pixel_scale, slice_scale, waves.min(), waves.max(), IPC, offset_par, offset_perp)):
        return PSFs[(pixel_scale, slice_scale, waves.min(), waves.max(), IPC, offset_par, offset_perp)], PSFs


    oneD_PSFs = zeros([len(waves), psfsize**2], dtype=float64)


    for j, wave in enumerate(waves):
        twoD_PSF = get_pixelized_PSF_noIPC(PSFs, pixel_scale = pixel_scale, slice_scale = slice_scale, wave = wave,
                                           offset_par = offset_par, offset_perp = offset_perp, psfsize = psfsize)
        
        for i in range(psfsize):
            oneD_PSFs[j, i*psfsize:(i+1)*psfsize] = twoD_PSF[:,i]

        #save_img(twoD_PSF, "2D.fits")
        #save_img(sliced, "sliced.fits")

    #save_img(oneD_PSFs, "oneD_PSFs.fits")
    oneD_PSFs_IPC = oneD_PSFs*(1 - 2.*IPC)
    oneD_PSFs_IPC[:,:-1] += oneD_PSFs[:,1:]*IPC
    oneD_PSFs_IPC[:,1:] += oneD_PSFs[:,:-1]*IPC
    #save_img(oneD_PSFs_IPC, "oneD_PSFs_IPC.fits")

    PSFs[(pixel_scale, slice_scale, waves.min(), waves.max(), IPC, offset_par, offset_perp)] = oneD_PSFs_IPC
    return oneD_PSFs_IPC, PSFs

"""

PSFs = initialize_PSFs(pixel_scales = [15], slice_scales = [30], PSF_source = "WebbPSF")

print PSFs
save_img(PSFs[(8126.0, 15, 30)], "tmp.fits")
aPSF = get_pixelized_PSF_noIPC(PSFs, pixel_scale = 15, slice_scale = 30, wave = 10000., offset_par = 0, offset_perp = 0)
save_img(aPSF, "tmp2.fits")
oneDPSF, NA = get_1D_pixelized_sliced_PSFs(PSFs, pixel_scale = 15, slice_scale = 30, waves = arange(5000., 20001., 100.), IPC = 0.02, offset_par = 0, offset_perp = 0)
save_img(oneDPSF, "tmp3.fits")

stop_here
"""

def image_to_cube(the_im, psfsize):
    return reshape(the_im, (len(the_im), psfsize, psfsize))

######################################### Signal-to-Noise Spectral Calculation #########################################


def photons_per_wave_to_flamb(photons_per_wave, meters2, waves, dwaves):
    erg_per_cm2 = 5.03411701e21 # meters^-3
    return photons_per_wave/(meters2*erg_per_cm2*(waves/1.e10) * dwaves)
    
def flamb_to_photons_per_wave(flamb, meters2, waves, dwaves):
    erg_per_cm2 = 5.03411701e21 # meters^-3
    return meters2*flamb*erg_per_cm2*(waves/1.e10) * dwaves


def get_thermal_background_per_pix_per_sec(waves, dwaves, pixel_scale, TTel, R = 1.2, throughput = 1., emiss = 0.05):
    """DR note 2/23/16: adding 1.36 fudge factor to match Chris Hirata. Not present in comparisons with Klaus and David.
    DR note 3/4/16: didn't run thermal to red enough wavelength; fudge factor 1.08."""

    return 1.08 * 4.43e28 * pixel_scale**2.* R**2. *dwaves * throughput * emiss / (
        waves**4. * (exp(1.43878e8/(TTel*waves)) - 1)
    )
    

#print get_thermal_background_per_pix_per_sec(waves = 19932.0001929, dwaves = 1., pixel_scale = 1., TTel = 282.)


"""
waves = arange(3000., 25000., 10.)
photons_per_arcsec_per_sec = get_thermal_background_per_pix_per_sec(waves = waves, dwaves = 1., pixel_scale = 1., TTel = 282.)
savetxt("output/thermal_background.txt", transpose(array([waves, photons_per_arcsec_per_sec,
                                                          photons_per_wave_to_flamb(photons_per_arcsec_per_sec, meters2, waves = waves, dwaves = 1.)])), header = "wavelength(A)\tphotons_per_arcsec2_sec_A", fmt = ["%f", "%g"])
flkdajfld
"""



#waves, dwaves = resolution_to_wavelengths(source_dir = "input", IFURfl = "IFU_R_160121.txt", min_wave = 4200., max_wave = 20000.)
#print len(waves)



def get_sncosmo(mdl, redshift, waves, fine_waves, phase, absmag = - 19.08):
    sncosmo_model = sncosmo.Model(source=mdl)
    sntype = {"s11-2005hl": "Ib", "s11-2005hm": "Ib", "s11-2006fo": "Ic", "s11-2006jo":
              "Ib", "hsiao": "Ia", "nugent-sn1a": "Ia", "salt2-extended": "Ia"}[mdl]
    cosmo = cosmology.FlatLambdaCDM(Om0 = 0.3, H0 = 70.)
    ten_pc_z = 2.33494867e-9
    assert abs(cosmo.distmod(z=ten_pc_z).value) < 1.e-3, "Distance modulus zeropoint wrong!"

    ampl = 1.
    for i in range(2):
        if not mdl.count("salt2"):
            sncosmo_model.set(z=ten_pc_z, t0=0., amplitude=ampl)
        else:
            sncosmo_model.set(z=ten_pc_z, t0=0., x0=ampl)

        mag = sncosmo_model.bandmag('bessellv', 'ab', 0.)
        print "mag, ampl ", mag, ampl
        ampl *= 10.**(0.4*(mag - absmag))
        

    if not mdl.count("salt2"):
        sncosmo_model.set(z=redshift, t0=0., amplitude=ampl)
    else:
        sncosmo_model.set(z=redshift, t0=0., x0=ampl)


    mu = cosmo.distmod(redshift).value
    print "mu ", mu

    f_lamb_SN = sncosmo_model.flux(phase*(1. + redshift), waves)*10.**(-0.4*mu)
    f_lamb_SN_fine = sncosmo_model.flux(phase*(1. + redshift), fine_waves)*10.**(-0.4*mu)
    return f_lamb_SN, f_lamb_SN_fine, sntype, sncosmo_model

def get_R07_noise(electron_count_rate, t_int, nframe, frame_time = 10.63, read_noise_per_frame = 10., mframe = 4.):
    group_time = t_int/(nframe - 1.)

    coeff1 = 12.*(nframe - 1) / (mframe*nframe*(nframe + 1))
    coeff2 = 1.2*(nframe**2 + 1)*(nframe - 1)/(nframe*(nframe + 1))
    coeff3 = -2*(mframe**2 - 1)*(nframe - 1)/(mframe*nframe*(nframe + 1))

    return sqrt(coeff1*read_noise_per_frame**2 + coeff2*group_time*electron_count_rate + coeff3*frame_time*electron_count_rate)


def get_spec_with_err(redshift, exp_time, phase = 0, gal_flamb = lambda x:0., pixel_scale = 0.075, slice_scale = 0.15, show_plots = 0,
                      dark_current = 0.01, offset_par = 0, offset_perp = 0,
                      mdl = 'hsiao', PSFs = None, photon_count_rates = None, output_dir = "output", othername = "",
                      IFURfl = "IFU_R_Content.txt", min_wave = 6000., max_wave = 20000.,
                      zodifl = "aldering.txt", effareafl = "IFU_effective_area_160513.txt", thermalfl = None,
                      source_dir = "input/", read_noise = None, white_noise = 15., read_noise_floor = 4., aper_rad_fn = None, waves = None, fine_waves = None,
                      restframe_bins = [3000., 4000., 5000., 6000., 8000., 10000., 12000.],
                      bad_pixel_rate = 0,
                      obsframe_bins = [7000, 8000., 10200, 12850, 16050, 20000.], TTel = 282., IPC = 0.02, nframe = None, use_R07_noise = False, psfsize = 7):

    if hasattr(effareafl, '__call__'):
        effective_meters2 = effareafl
    else:
        effective_meters2 = interpfile(source_dir + "/" + effareafl, bounds_error = False)

    if hasattr(zodifl, '__call__'):
        zodi_flamb = zodifl
    else:
        zodi_flamb = interpfile(source_dir + "/" + zodifl)
    #galaxy_flamb = interpfile(source_dir + "/" + "NGC_7591_beta_spec.dat", norm = True)

   
    AB_mags = [0., 20., 25.] # AB mags to generate

    
    
    if read_noise == None:
        read_noise = (read_noise_floor, sqrt(
            12.* white_noise**2. * (exp_time/0.65 - 1.)/ (exp_time/0.65) / (exp_time/0.65  + 1.)
        ))

    # Starting model:
    
    waves, dwaves = resolution_to_wavelengths(source_dir, IFURfl, min_wave, max_wave, waves)

    print "waves ", waves.min(), waves.max(), len(waves)
    if show_plots:
        savetxt("output/wavelengths.txt", zip(arange(1, len(waves + 1)), waves), fmt = ["%i", "%f"])
    if fine_waves == None:
        fine_waves = arange(3000., 22001., 10.)

    try:
        mdl_len = len(mdl)
    except:
        mdl_len = 0

        
    if mdl_len == len(waves):
        f_lamb_SN = mdl
        f_lamb_SN_fine = interp1d(waves, mdl, bounds_error = False, fill_value = 0)(fine_waves)
        sntype = "Unknown"
        mdl = "Unknown"
    elif glob.glob(source_dir + "/" + mdl) == []:
        f_lamb_SN, f_lamb_SN_fine, sntype, NA = get_sncosmo(mdl, redshift, waves, fine_waves, phase)

    else:
        # Okay. Reading from file.

        flambfn = interpfile(source_dir + "/" + mdl)
        mu = None
        sntype = "NA"
        
        f_lamb_SN = flambfn(waves)
        f_lamb_SN_fine = flambfn(fine_waves)
        f_lamb_SN_at_max = flambfn(waves)

    if PSFs == None:
        PSFs = initialize_PSFs(pixel_scales = [int(round(pixel_scale/0.005))], slice_scales = [int(round(slice_scale/0.005))])

    # For reference: waves, dwaves = resolution_to_wavelengths(source_dir, IFURfl, min_wave, max_wave, waves)

    oneD_PSFs, PSFs = get_1D_pixelized_sliced_PSFs(PSFs, pixel_scale = int(round(pixel_scale/0.005)), slice_scale = int(round(slice_scale/0.005)), waves = waves, IPC = IPC, offset_par = offset_par, offset_perp = offset_perp, psfsize=psfsize)
    
    #print "f_lamb_SN ", f_lamb_SN

    try:
        plt_root = "%s_z=%.2f_exp=%.1f_ph=%.1f_gal=%.1g_PSF=%s_pxsl=%.3f_ssl=%i_mdl=%s_type=%s_%s_zodi=%s_opar=%iperp=%i_RN=%.1f" % (othername, redshift, exp_time, phase, gal_flamb(15000.), PSFs["PSF_source"],
                                                                                                                                     pixel_scale, slice_scale, mdl.split("/")[-1].split(".")[0], sntype, IFURfl.split("/")[-1].split(".")[0], zodifl.split("/")[-1].split(".")[0],
                                                                                                                                     offset_par, offset_perp, sqrt(read_noise[0]**2 + read_noise[1]**2))
    except:
        plt_root = ""

    if show_plots:
        savetxt(output_dir + "/SN_template_z=%.2f_ph=%.1f_mdl=%s_type=%s.txt" % (redshift, 0., mdl.split("/")[-1].split(".")[0], sntype),
                transpose(array(  [fine_waves, f_lamb_SN_fine]   )), fmt = ["%f", "%g"])
    photons_SN_per_sec = flamb_to_photons_per_wave(f_lamb_SN, effective_meters2(waves), waves, dwaves)

    if photon_count_rates != None:
        # This is a dictionary of wavelength range, photon_count_rate
        for waverange in photon_count_rates:
            inds = where((waves >= waverange[0])*(waves < waverange[1]))
            norm_factor = photon_count_rates[waverange]/sum(photons_SN_per_sec[inds])
            print waverange, "norm_factor", norm_factor
            f_lamb_SN[inds] *= norm_factor
            f_lamb_SN_at_max[inds] *= norm_factor
            photons_SN_per_sec[inds] *= norm_factor

        model = None # Prevent the now incorrect model from being used again

    zodi_photons_per_sec_perarcsec2 = flamb_to_photons_per_wave(10.**(zodi_flamb(waves)), effective_meters2(waves), waves, dwaves)
    galaxy_photons_per_sec_perarcsec2 = flamb_to_photons_per_wave(gal_flamb(waves), effective_meters2(waves), waves, dwaves)


    if thermalfl == None:
        thermal_background_per_pix_per_sec = get_thermal_background_per_pix_per_sec(waves, dwaves, sqrt(pixel_scale*slice_scale), TTel = TTel,
                                                                                    throughput = effective_meters2(waves)/(pi*1.2**2.))
        
        thermal_flamb_arcsec2 = photons_per_wave_to_flamb(get_thermal_background_per_pix_per_sec(waves = waves, dwaves = 1., pixel_scale = 1.,
                                                                                                 TTel = TTel, throughput = effective_meters2(waves)/(pi*1.2**2.)),
                                                          effective_meters2(waves), waves = waves, dwaves = 1.) # For pixel scale 1, get thermal background
    else:
        thermal_flamb_arcsec2 = interpfile(source_dir + "/" + thermalfl)
        thermal_flamb_arcsec2 = 10.**(thermal_flamb_arcsec2(waves))
        thermal_background_per_pix_per_sec = flamb_to_photons_per_wave(thermal_flamb_arcsec2, effective_meters2(waves), waves, dwaves)*pixel_scale*slice_scale


    AB_mag_per_sec = {}
    for AB_mag in AB_mags:
        AB_mag_per_sec[AB_mag] = flamb_to_photons_per_wave(0.10884806248 * 10**(-0.4*AB_mag) /waves**2., effective_meters2(waves), waves, dwaves)


    
    SN_photon_image = transpose(transpose(oneD_PSFs)*photons_SN_per_sec)*exp_time
    zodi_image = transpose(transpose(zeros(oneD_PSFs.shape, dtype=float64)) + zodi_photons_per_sec_perarcsec2*exp_time*pixel_scale*slice_scale)
    gal_image = transpose(transpose(zeros(oneD_PSFs.shape, dtype=float64)) + galaxy_photons_per_sec_perarcsec2*exp_time*pixel_scale*slice_scale)
    thermal_image = transpose(transpose(zeros(oneD_PSFs.shape, dtype=float64)) + thermal_background_per_pix_per_sec*exp_time)
    dark_current_image = zeros(oneD_PSFs.shape, dtype=float64) + dark_current*exp_time
    #dark_current_image = random.lognormal(log(dark_current), 1, size = oneD_PSFs.shape)*exp_time
    


    AB_mag_images = {}
    for AB_mag in AB_mags:
        AB_mag_images[AB_mag] = transpose(transpose(oneD_PSFs)*AB_mag_per_sec[AB_mag])*exp_time

    total_image = SN_photon_image + zodi_image + gal_image + thermal_image + dark_current_image
    #assert all(total_image < 6.e4), "Saturated pixels found!"

    if not use_R07_noise:
        total_noise = sqrt(total_image + read_noise[0]**2. + read_noise[1]**2.)
    else:
        total_noise = get_R07_noise(electron_count_rate = total_image/exp_time, t_int = exp_time, nframe = nframe)

    if bad_pixel_rate > 0:
        total_noise += (random.random(size = total_noise.shape) < bad_pixel_rate)*total_noise[0,0]*100.


    sim_noise = random.normal(size = total_noise.shape)*total_noise
    sim_weight = 1./total_noise**2.


    SN_image_with_noise = SN_photon_image + sim_noise

 
    SN_with_noise_cube = image_to_cube(SN_image_with_noise, psfsize)
    SN_with_noise_mean2d = mean(SN_with_noise_cube, axis = 0)

    """
    # Make the cubes. Signal cubes should be scaled by 1/slice_in_pixels to preserve flux. Noise cubes should be scaled by 1/sqrt(sice_in_pixels) to preserve noise.
    SN_with_noise_cube = image_to_cube(SN_image_with_noise, slice_in_pixels)
    total_noise_cube = image_to_cube(total_noise, slice_in_pixels)
    SN_photon_cube = image_to_cube(SN_photon_image, slice_in_pixels)
    PSF_cube = image_to_cube(oneD_PSFs, slice_in_pixels)"""
    
    #SN_photon_mean2d = mean(SN_photon_cube, axis = 0)

    """
    sn_sep = sep.extract(SN_photon_mean2d, thresh = 0.001)
    SN_aperture_StoN = []
    SN_aperture_signal = []
    aper_rad_by_wave = []
    aper_EE = []

    for i, wave in enumerate(waves):
        PSFsig = PSFs[(int(round(pixel_scale/0.005)), "sigmaFN")](wave)
        if aper_rad_fn == None:
            aper_rad = 2.5*PSFsig/float(pixel_scale/0.005)
        else:
            aper_rad = aper_rad_fn(wave)/float(pixel_scale)
            
        signal, noise, NA = sep.sum_circle(SN_photon_cube[i], x = sn_sep["x"], y = sn_sep["y"], r = aper_rad, err = total_noise_cube[i], subpix = 0)
        apsum, NA, NA = sep.sum_circle(PSF_cube[i], x = sn_sep["x"], y = sn_sep["y"], r = aper_rad, subpix = 0)
        
        try:
            aper_EE.append(apsum[0])
            SN_aperture_StoN.append(signal[0]/noise[0])
            SN_aperture_signal.append(signal[0])
        except:
            aper_EE.append(-1)
            SN_aperture_StoN.append(-1)
            SN_aperture_signal.append(-1)
        aper_rad_by_wave.append(aper_rad)
    """

    if show_plots:
        save_img(oneD_PSFs, output_dir + "/oneD_PSFs_" + plt_root + ".fits")
        save_img(SN_photon_image, output_dir + "/SN_image_" + plt_root + ".fits")
        save_img(zodi_image, output_dir + "/zodi_image_" + plt_root + ".fits")
        save_img(gal_image, output_dir + "/gal_image_" + plt_root + ".fits")
        save_img(thermal_image, output_dir + "/thermal_image_" + plt_root + ".fits")
        save_img(total_image, output_dir + "/total_image_" + plt_root + ".fits")
        save_img(total_noise, output_dir + "/total_noise_" + plt_root + ".fits")
        save_img(total_image + sim_noise, output_dir + "/image_with_noise_" + plt_root + ".fits")
        save_img(SN_with_noise_cube, output_dir + "/SN_with_noise_" + plt_root + "_cube.fits")
        save_img(SN_with_noise_mean2d, output_dir + "/SN_with_noise_" + plt_root + "_mean2d.fits")

        plt.plot(waves, PSFs[(int(round(pixel_scale/0.005)), int(round(slice_scale/0.005)), "sigmaFN")](waves))
        plt.savefig(output_dir + "/psfsigma_" + plt_root + ".pdf", bbox_inches = 'tight')
        plt.close()
    

    extract_denom = sum(sim_weight*oneD_PSFs*oneD_PSFs, axis = 1)
    PSF_enclosed_flux = sum(oneD_PSFs, axis = 1)

    extracted_SN = sum(sim_weight*SN_image_with_noise*oneD_PSFs, axis = 1)/extract_denom
    extracted_noise = 1./sqrt(extract_denom)
    extracted_zodinoise = sqrt(sum(sim_weight*zodi_image*oneD_PSFs, axis = 1)/extract_denom/PSF_enclosed_flux)
    extracted_darknoise = sqrt(sum(sim_weight*dark_current_image*oneD_PSFs, axis = 1)/extract_denom/PSF_enclosed_flux)
    extracted_readfloor = sqrt(sum(sim_weight*read_noise[0]**2. * oneD_PSFs, axis = 1)/extract_denom/PSF_enclosed_flux)
    extracted_readwhite = sqrt(sum(sim_weight*read_noise[1]**2. * oneD_PSFs, axis = 1)/extract_denom/PSF_enclosed_flux)
    extracted_thermalnoise = sqrt(sum(sim_weight*thermal_image*oneD_PSFs, axis = 1)/extract_denom/PSF_enclosed_flux)
    extracted_galnoise = sqrt(sum(sim_weight*gal_image*oneD_PSFs, axis = 1)/extract_denom/PSF_enclosed_flux)
    extracted_SNnoise = sqrt(sum(sim_weight*SN_photon_image*oneD_PSFs, axis = 1)/extract_denom/PSF_enclosed_flux)

    eff_noise = (extracted_readfloor**2. + extracted_readwhite**2.)/(read_noise[0]**2. + read_noise[1]**2.)



    yerr_flamb = photons_per_wave_to_flamb(extracted_noise/exp_time, effective_meters2(waves), waves, dwaves)
    yvalswitherr_flamb = photons_per_wave_to_flamb(extracted_SN/exp_time, effective_meters2(waves), waves, dwaves)
    
    
    signal_to_noises = {"obs_frame_band": {}, # band-integrated
                        "obs_frame": {}, # per resolution element
                        "rest_frame_band_S/N": {},
                        "rest_frame_mean_S/N": {},
                        "thermal_flamb_arcsec-2": thermal_flamb_arcsec2,
                        "plt_root": plt_root,
                        "PSFs": PSFs,
                        "PSF_enclosed_flux": PSF_enclosed_flux,
                        "rest_frame_band_mag": {},
                        "photons/s_obs": {}, "dark_photnoise": {}, "n_reselmnts": {}, "thermal_photnoise": {}, "zodi_photnoise": {}, "read_noise_white": {}, "read_noise_floor": {}, "SN_photnoise": {},
                        "f_lamb_SN": f_lamb_SN, "extracted_SN": extracted_SN,
                        "obs_waves": waves, "obs_dwaves": dwaves,
                        "spec_S/N": f_lamb_SN/yerr_flamb,# "spec_aperture_StoN": array(SN_aperture_StoN), "spec_aperture_signal": array(SN_aperture_signal),
                        #"aper_rad_by_wave": aper_rad_by_wave, "aperture_EE": aper_EE,
                        "PSF_wghtd_e_SN": sum(oneD_PSFs*SN_photon_image, axis = 1),
                        "PSF_wghtd_e_allbutdark": sum(oneD_PSFs*(SN_photon_image + zodi_image + gal_image + thermal_image), axis = 1),
                        "phot/s_image": (SN_photon_image + zodi_image + gal_image + thermal_image)/exp_time,
                        "noise_sources": {"total": extracted_noise, "zodi": extracted_zodinoise, "thermal": extracted_thermalnoise,
                                          "dark": extracted_darknoise, "readwhite": extracted_readwhite, "readfloor": extracted_readfloor, "SN": extracted_SNnoise, "galaxy": extracted_galnoise}
                        }



    for i in range(len(obsframe_bins) - 1):
        inds = where((waves >= obsframe_bins[i])*(waves < obsframe_bins[i+1]))
        obskey = (obsframe_bins[i], obsframe_bins[i+1])

        signal_to_noises["n_reselmnts"][obskey] = float(len(inds[0]))*0.5
        signal_to_noises["dark_photnoise"][obskey] = sqrt(dot(extracted_darknoise[inds], extracted_darknoise[inds])) # Quadrature sum
        signal_to_noises["read_noise_white"][obskey] = sqrt(dot(extracted_readwhite[inds], extracted_readwhite[inds])) # Quadrature sum
        signal_to_noises["read_noise_floor"][obskey] = sqrt(dot(extracted_readfloor[inds], extracted_readfloor[inds])) # Quadrature sum
        signal_to_noises["zodi_photnoise"][obskey] = sqrt(dot(extracted_zodinoise[inds], extracted_zodinoise[inds])) # Quadrature sum
        signal_to_noises["thermal_photnoise"][obskey] = sqrt(dot(extracted_thermalnoise[inds], extracted_thermalnoise[inds])) # Quadrature sum
        signal_to_noises["SN_photnoise"][obskey] = sqrt(dot(extracted_SNnoise[inds], extracted_SNnoise[inds])) # Quadrature sum

        this_errs = sqrt(sum((yerr_flamb[inds]*dwaves[inds]*waves[inds])**2.))
        this_photons = sum(f_lamb_SN[inds]*dwaves[inds]*waves[inds])
        

        #weighted_photons = sum(these_weights*these_photons)/sum(these_weights)
        #weighted_errs = 1./sqrt(sum(these_weights))

        nresl = 0.5*float(len(inds[0]))

        print "MeanS/N", obsframe_bins[i], obsframe_bins[i+1], this_photons/this_errs, this_photons/this_errs/sqrt(nresl)
        signal_to_noises["obs_frame"][obskey] = this_photons/this_errs/sqrt(nresl)
        signal_to_noises["obs_frame_band"][obskey] = this_photons/this_errs
        signal_to_noises["photons/s_obs"][obskey] = sum(photons_SN_per_sec[inds])


    




    s_to_n_tmp = integrate_spec(waves, dwaves, redshift, restframe_bins, f_lamb_SN, yerr_flamb)
    
    for key in s_to_n_tmp:
        signal_to_noises[key] = s_to_n_tmp[key]            


    if show_plots:

        binwaves = arange(min(waves), max(waves), 150.) #clip(exp(arange(log(min(waves)), log(max(waves)), 1./80.)), min(waves), max(waves))
        f_lamb_SN_bin, f_lamb_SN_binerr = bin_vals_fixed_bins(waves, yvalswitherr_flamb, yerr_flamb, binwaves)
        plt.plot(waves, f_lamb_SN, color = 'k', zorder = 1)
        plt.errorbar(binwaves, f_lamb_SN_bin, yerr = f_lamb_SN_binerr, fmt = '.', label = plt_root, capsize = 0,
                     color = 'b'*(1 - dark_background) + dark_background*'w')
        plt.legend(loc = 'best', fontsize = 8)
        plt.ylabel("$f_{\lambda}$")
        plt.xlabel("Observer-Frame Wavelength (Binned)")
        ylim = list(plt.ylim())
        ylim[0] = 0
        plt.ylim(ylim)
        plt.savefig(output_dir + "/spec_binned_" + plt_root + ".pdf", bbox_inches = 'tight')
        plt.close()


        plt.plot(waves, eff_noise)
        plt.savefig(output_dir + "/eff_noise_" + plt_root + ".pdf")
        plt.close()


        plt.errorbar(waves, yvalswitherr_flamb, yerr = yerr_flamb, fmt = '.', label = plt_root, capsize = 0)
        plt.legend(loc = 'best', fontsize = 8)
        plt.ylabel("$f_{\lambda}$")
        plt.xlabel("Observer-Frame Wavelength")
        plt.ylim(ylim)
        plt.savefig(output_dir + "/spec_" + plt_root + ".pdf", bbox_inches = 'tight')
        plt.close()

        plt.plot(waves, f_lamb_SN/yerr_flamb, '.', color = 'k'*(1 - dark_background) + 'w'*dark_background, label = plt_root)
        plt.ylabel("Signal-to-Noise per Wavelength (half resolution element)")
        plt.xlabel("Observer-Frame Wavelength")
        plt.legend(loc = 'best', fontsize = 8)
        plt.savefig(output_dir + "/spec_StoN_" + plt_root + ".pdf", bbox_inches = 'tight')
        plt.close()


        

        #plt.plot(waves, signal_to_noises["spec_aperture_StoN"], '.', color = 'k'*(1 - dark_background) + 'w'*dark_background, label = plt_root)
        #plt.ylabel("Aperture Signal-to-Noise per Wavelength (half resolution element)")
        #plt.xlabel("Observer-Frame Wavelength")
        #plt.legend(loc = 'best', fontsize = 8)
        #plt.savefig(output_dir + "/spec_aperture_StoN_" + plt_root + ".pdf", bbox_inches = 'tight')
        #plt.close()
        

        savetxt(output_dir + "/sim_spectrum_" + plt_root + ".txt", transpose(array([waves, f_lamb_SN, yerr_flamb])), header = "#waves   template   err", fmt = ["%f", "%g", "%g"])
        savetxt(output_dir + "/sim_noise_" + plt_root + ".txt",  transpose(array([waves] + 
                 [photons_per_wave_to_flamb(item/exp_time, effective_meters2(waves), waves, dwaves) for item in [extracted_zodinoise, extracted_darknoise, extracted_readfloor, extracted_readwhite, extracted_SNnoise]]
                                                                             )), header = "#waves   zodi   dark   readfloor  readwhite   SN", fmt = ["%f"] + ["%g"]*5)

    return signal_to_noises


def solve_for_exptime(S_to_N, redshift, PSFs, key1 = "obs_frame", key2 = (10200, 12850), **kwargs):
    """Solves for exposure times. Example kwargs include: pixel_scale = 0.075, slice_in_pixels = 2, dark_current = 0.01, photon_count_rates = None, gal_norm = 0, mdl = 'hsiao'"""

    times = [1000.]
    SNs = [get_spec_with_err(redshift, exp_time = item, PSFs = PSFs, show_plots = 0, **kwargs)[key1][key2]
           for item in times]
    
    if SNs[0] == None:
        return None

    while (max(SNs) < S_to_N) or (min(SNs) > S_to_N):
        if max(SNs) < S_to_N:
            times.append(max(times)*4.)
        if min(SNs) > S_to_N:
            times.append(min(times)/4.)
        SNs.append(get_spec_with_err(redshift, exp_time = times[-1], PSFs = PSFs, show_plots = 0, **kwargs)[key1][key2]
        )

    last_guess = 0.2
    guess = 0.1
    while abs(log(last_guess/guess)) > 0.001:
        # While guesses are more than 0.1% apart:
        last_guess = guess
        guess = interp1d(SNs, times)(S_to_N)
        print "guess ", guess
        times.append(guess)
        SNs.append(get_spec_with_err(redshift, exp_time = times[-1], PSFs = PSFs, show_plots = 0, **kwargs)[key1][key2]
        )
    guess = interp1d(SNs, times)(S_to_N)
    print "Returning ", guess
    
    return guess


######################################### Signal-to-Noise Imaging Calculation #########################################


def get_imaging_SN(PSFs, exp_time, effective_meters2_fl, wavemin = 4000, wavemax = 25000, waves = None, redshift=0, phase=0, gal_flamb = lambda x:0., pixel_scale = 0.11, IPC = 0.02, offset_par = 5, offset_perp = 5, source_dir = "input", zodi_fl = "aldering.txt", mdl = "hsiao", dark_current = 0.015, TTel = 282., verbose = False, approximate_PSF = True, bad_pixel_rate = 0, read_noise_floor = 5., read_noise_white = 20.):
    scale = int(round(pixel_scale/0.005))
    if any(waves == None):
        waves = arange(wavemin, wavemax, 50.) # These should span (and perhaps slightly overfill) the filter
        dwaves = ones(len(waves), dtype=float64)*50.
    else:
        dwaves = waves[1:] - waves[:-1]
        dwaves = concatenate(([dwaves[0]], dwaves))

    if hasattr(effective_meters2_fl, '__call__'):
        effective_meters2 = effective_meters2_fl
    else:
        effective_meters2 = interpfile(source_dir + "/" + effective_meters2_fl, bounds_error = False)

    if hasattr(zodi_fl, '__call__'):
        zodi_flamb = zodi_fl
    else:
        zodi_flamb = interpfile(source_dir + "/" + zodi_fl)


    
    """    
    read_noise = sqrt(5.**2 + 
                      3.* 20.**2. / (exp_time/11.3)
                 )
    """
    read_noise = sqrt(read_noise_floor**2 + 
                      12.* read_noise_white**2. * (exp_time/2.825 - 1.)/ (exp_time/2.825) / (exp_time/2.825  + 1.)
                 )


    if verbose:
        print "read_noise", read_noise

    try:
        model_len = len(mdl)
    except:
        model_len = None

    if mdl == 'hsiao':
        cosmo = cosmology.FlatLambdaCDM(Om0 = 0.3, H0 = 70.)
        ten_pc_z = 2.33494867e-9
        assert abs(cosmo.distmod(z=ten_pc_z).value) < 1.e-3, "Distance modulus zeropoint wrong!"


        model = sncosmo.Model(source=mdl)
        sntype = {"s11-2005hl": "Ib", "s11-2005hm": "Ib", "s11-2006fo": "Ic", "s11-2006jo": "Ib", "hsiao": "Ia", "nugent-sn1a": "Ia", "salt2-extended": "Ia"}[mdl]
        

        ampl = 1.
        for i in range(3):
            model.set(z=ten_pc_z, t0=0., amplitude=ampl)
            
            mag = model.bandmag('bessellv', 'ab', 0.)
            if verbose:
                print "mag, ampl ", mag, ampl
            ampl *= 10.**(0.4*(mag -- 19.08))

        model.set(z=redshift, t0=0., amplitude=ampl)

        mu = cosmo.distmod(redshift).value
        if verbose:
            print "mu ", mu
        
        f_lamb_SN = model.flux(phase*(1. + redshift), waves)*10.**(-0.4*mu)
        f_lamb_SN_at_max = model.flux(0., waves)*10.**(-0.4*mu)
    elif type(mdl) == types.FunctionType:
        f_lamb_SN = mdl(waves)
        f_lamb_SN_at_max = mdl(waves)
        sntype = "None"
    elif model_len == len(waves):
        f_lamb_SN = mdl
        f_lamb_SN_at_max = mdl
        sntype = "None"
    elif len(glob.glob(source_dir + "/" + mdl)) == 1:
        mdlfn = interpfile(source_dir + "/" + mdl)
        f_lamb_SN = mdlfn(waves)
        f_lamb_SN_at_max = mdlfn(waves)
        sntype = "None"

    if 0:
        plt_root = "z=%.2f_exp=%.1f_ph=%.1f_gal=%.1g_PSF=%s_pxsl=%.3f_mdl=%s_type=%s_%s" % (redshift, exp_time, phase, gal_flamb(15000.), PSFs["PSF_source"], pixel_scale, mdl, sntype, effective_meters2_fl)

    #if show_plots:
    #    writecol(output_dir + "/SN_template_" + plt_root + ".txt", transpose(array([arange(3000., 25001., 10.), model.flux(0., arange(3000., 25001., 10.))*10.**(-0.4*mu)])))
    photons_SN_per_secwave = flamb_to_photons_per_wave(f_lamb_SN, effective_meters2(waves), waves, dwaves)

    zodi_photons_per_pixsec = sum(flamb_to_photons_per_wave(10.**(zodi_flamb(waves)), effective_meters2(waves), waves, dwaves))*pixel_scale**2.
    if verbose:
        print "zodi_photons_per_pixsec ", effective_meters2_fl, zodi_photons_per_pixsec
    galaxy_photons_per_pixsec = sum(flamb_to_photons_per_wave(gal_flamb(waves), effective_meters2(waves), waves, dwaves))*pixel_scale**2.
    thermal_photons_per_pixsec = sum(get_thermal_background_per_pix_per_sec(waves, dwaves, pixel_scale = pixel_scale,
                                                                            TTel = TTel, throughput = effective_meters2(waves)/(pi*1.2**2.))
                                     )

    if verbose:
        print "thermal_photons_per_pixsec ", effective_meters2_fl, thermal_photons_per_pixsec

    if PSFs == None:
        PSFs = initialize_PSFs(pixel_scales = [int(round(pixel_scale/0.005))], slice_scales = [int(round(pixel_scale/0.005))])
    if approximate_PSF:
        the_PSF = get_approximate_pixelized_broadband_PSF(PSFs, scale = int(round(pixel_scale/0.005)), waves = waves, weights = photons_SN_per_secwave/sum(photons_SN_per_secwave), IPC = IPC, offset_par = offset_par, offset_perp = offset_perp)
    else:
        the_PSF = get_pixelized_broadband_PSF(PSFs, scale = int(round(pixel_scale/0.005)), waves = waves, weights = photons_SN_per_secwave/sum(photons_SN_per_secwave), IPC = IPC, offset_par = offset_par, offset_perp = offset_perp)
        

    SN_photon_image = the_PSF*sum(photons_SN_per_secwave)*exp_time
    if verbose:
        print "SN photons/sec ", redshift, sum(photons_SN_per_secwave)
    zodi_image = ones(the_PSF.shape, dtype=float64)*zodi_photons_per_pixsec*exp_time
    gal_image = ones(the_PSF.shape, dtype=float64)*galaxy_photons_per_pixsec*exp_time
    thermal_image = ones(the_PSF.shape, dtype=float64)*thermal_photons_per_pixsec*exp_time
    dark_current_image = ones(the_PSF.shape, dtype=float64)*dark_current*exp_time

    total_image = SN_photon_image + zodi_image + gal_image + thermal_image + dark_current_image
    total_noise = sqrt(total_image + read_noise**2.)

    if bad_pixel_rate > 0:
        total_noise += (random.random(size = total_noise.shape) < bad_pixel_rate)*total_noise[0,0]*100.

    sim_weight = 1./total_noise**2.

    #print SN_photon_image
    #print the_PSF
    #print sim_weight
    #print "exp_time ", exp_time



    signal_to_noise = sum(SN_photon_image*the_PSF*sim_weight)/sqrt(sum(the_PSF**2. * sim_weight))
    ETC_results = {}
    ETC_results["PSF_phot_S/N"] = signal_to_noise
    ETC_results["SN_photons/s"] = sum(photons_SN_per_secwave)
    AB_0 = sum(flamb_to_photons_per_wave(0.10884806248/waves**2., effective_meters2(waves), waves, dwaves))
    ETC_results["AB_mag"] = -2.5*log10(ETC_results["SN_photons/s"]/AB_0)
    ETC_results["zodi/s"] = zodi_photons_per_pixsec
    ETC_results["thermal/s"] = thermal_photons_per_pixsec
    

    return ETC_results

    
def emulate_ETC(**kwargs):
    sample_result = get_spec_with_err(**kwargs)
    


#get_spec_with_err(redshift = 0., exp_time = 900., show_plots = 1, phase = 0, mdl = "input/5000K_21st.txt", min_wave = 4000., othername = "tmp1_")
#get_spec_with_err(redshift = 0., exp_time = 900., show_plots = 1, phase = 0, mdl = "input/5000K_21st.txt", min_wave = 4000., othername = "tmp2_", zodifl = "thermal.txt")
#get_spec_with_err(redshift = 0., exp_time = 900., show_plots = 1, phase = 0, mdl = "input/5000K_21st.txt", min_wave = 4000., othername = "tmp3_", thermalfl = "thermal.txt")


if __name__ == "__main__":
    print "Running as a demo!"

    """
    print "Reading PSFs..."
    PSFs = initialize_PSFs(scales = [30,60], PSF_source = "exp") # Native is in units of 5 mas, so this is 0".075, 0".11, and 0".15
    
    for pixel_scale in [0.15, 0.3]:
        args = {"gal_flamb": lambda x:0., "pixel_scale": pixel_scale, "slice_in_pixels": 1, "dark_current": 0.01, "mdl": 'hsiao', "PSFs": PSFs, "IFURfl": "IFU_R_Content.txt", "min_wave": 10000.}
        PSFs = get_spec_with_err(redshift = 2.0, exp_time = 3600., show_plots = 1, phase = 0, **args)["PSFs"]
    stop_here
    """

    PSFs = initialize_PSFs(pixel_scales = [10, 15, 22, 30], slice_scales = [30, 30, 22, 30], PSF_source = "WebbPSF") # Native is in units of 5 mas, so this is 0".075, 0".11, and 0".15

    # 6e-19
    args = {"gal_flamb": lambda x:0, "pixel_scale": 0.05, "slice_scale": 0.15, "dark_current": 0.003,
            "mdl": 'hsiao', "PSFs": PSFs, "IFURfl": "IFU_R_160720.txt", "min_wave": 4200., "max_wave": 21000., "bad_pixel_rate": 0.0, "offset_par": 4}
    PSFs = get_spec_with_err(redshift = 1.0, exp_time = 1089.2, show_plots = 1, phase = 0, **args)["PSFs"]


    #ETC_result = get_spec_with_err(redshift = 1.0, exp_time = 1089., phase = 0, show_plots = 1, **args)
    #args["bad_pixel_rate"] = 0.01
    fffff

    for redshift in [0.5, 1.0, 1.5]:
        exptime = solve_for_exptime(S_to_N = 3.5*sqrt(15.), redshift = redshift, key1 = "rest_frame_band_S/N", key2 = (5000, 6000), phase = 0, **args)
        print "# redshift/ exposure time ", redshift, exptime

        for phase in [-8, -6, 6, 8, 10, 12]:
            ETC_result = get_spec_with_err(redshift = redshift, exp_time = exptime, phase = phase, show_plots = 0, **args)
            print "#  ", redshift, exptime, phase, ETC_result["rest_frame_band_S/N"][(5000, 6000)]


    stop_here_now

    """
    import cProfile, pstats, StringIO
    pr = cProfile.Profile()
    pr.enable()

    get_spec_with_err(redshift = 1.0, exp_time = 1000., show_plots = 0, phase = 0, **args)

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()
    fffff
    """

    #for redshift in arange(0.15, 1.66, 0.1):
    #    exptime = solve_for_exptime(S_to_N = 14.14, redshift = redshift, key1 = "rest_frame_mean_S/N", key2 = (3900, 4900), restframe_bins = [3900, 4900], phase = 0, **args)
    #    print "##redshift/ exposure time ", redshift, exptime

    """
    args["mdl"] = "27th_mag_ST.txt"
    args["psfsize"] = 11
    args["PSFs"] = None
    args["min_wave"] = 4200.
    args["max_wave"] = 21000.

    get_spec_with_err(redshift = 0., exp_time = 100000, show_plots = 1, phase = 0, **args)

    stop_here
    """

    for redshift in [0.05, 0.5, 1.0, 1.5]:
        exptime = solve_for_exptime(S_to_N = 10., redshift = redshift, key1 = "obs_frame", key2 = (10200, 12850), phase = 0, **args)
        print "redshift/ exposure time ", redshift, exptime

        get_spec_with_err(redshift = redshift, exp_time = exptime, show_plots = 1, phase = 0, **args)

        """
        for offset_perp in arange(15, 17.):
            exptime = solve_for_exptime(S_to_N = 10., redshift = redshift, key1 = "obs_frame", key2 = (10200, 12850), phase = 0, offset_perp = offset_perp, **args)
            get_spec_with_err(redshift = redshift, exp_time = exptime, show_plots = 1, phase = 0, offset_perp = offset_perp, **args)
        """
    stop_here


    if 1:
        #import cProfile, pstats, StringIO
        #pr = cProfile.Profile()
        #pr.enable()
        
        t = time.time()
        ETC_result = get_imaging_SN(PSFs, exp_time = 100., redshift = 1.0, effective_meters2_fl = "F184.txt", waves = arange(4000., 25000., 20.), approximate_PSF = False)
        print time.time() - t
        print ETC_result

        t = time.time()
        ETC_result = get_imaging_SN(PSFs, exp_time = 100., redshift = 1.0, effective_meters2_fl = "F184.txt", waves = arange(4000., 25000., 20.), approximate_PSF = True)
        print time.time() - t
        print ETC_result

        
        #pr.disable()
        #s = StringIO.StringIO()
        #sortby = 'cumulative'
        #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        #ps.print_stats()
        #print s.getvalue()
    

    get_spec_with_err(redshift = 0.1, exp_time = 30, show_plots = 1, phase = 0, **args)

    get_spec_with_err(redshift = 0., exp_time = 900., show_plots = 1, phase = 0, mdl = "5000K_21st.txt", min_wave = 4000.)
    get_spec_with_err(redshift = 0., exp_time = 9., show_plots = 1, phase = 0, mdl = "../../../../calspec/gd71_stisnic_006.ascii", min_wave = 4000.)

    get_spec_with_err(redshift = 0., exp_time = 1.3, show_plots = 1, phase = 0, mdl = "../../../../calspec/bd_17d4708_stisnic_006.ascii", min_wave = 4000.)


    args["zodifl"] = "zodi_medium.txt"
    exptime = solve_for_exptime(S_to_N = 10., redshift = 1.4, key1 = "rest_frame_band_S/N", key2 = (5000, 6000), phase = 0, **args)
    print "exp time", exptime

    args["zodifl"] = "aldering.txt"
    exptime = solve_for_exptime(S_to_N = 10., redshift = 1.4, key1 = "rest_frame_band_S/N", key2 = (5000, 6000), phase = 0, **args)
    print "exp time", exptime



    print "GRISM!!!!"
    args = {"gal_flamb": lambda x:0., "pixel_scale": 0.11, "slice_in_pixels": 1, "dark_current": 1.0, "mdl": 'hsiao', "PSFs": PSFs, "IFURfl": "IFU_R_Content.txt", "min_wave": 4000.}
    #PSFs = initialize_PSFs(scales = [15, 22, 30], PSF_source = "WebbPSF") # Native is in units of 5 mas, so this is 0".075, 0".11, and 0".15

    exptime = solve_for_exptime(S_to_N = 10., redshift = 1.0, key1 = "obs_frame", key2 = (10200, 12850), phase = 0, **args)
    get_spec_with_err(redshift = 1.0, exp_time = exptime, show_plots = 1, phase = 0, **args)
