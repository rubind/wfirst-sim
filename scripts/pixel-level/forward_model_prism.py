import numpy as np
from DavidsNM import save_img, miniLM_new
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
from pixel_level_ETC2 import resolution_to_wavelengths, initialize_PSFs, photons_per_wave_to_flamb, flamb_to_photons_per_wave
from FileRead import file_to_fn
import os
import sys
wfirst_path = os.environ["WFIRST"]
wfirst_data_path = os.environ["WFIRST_SIM_DATA"]




def pullfn_per_pix(P, all_data):
    ifn = interp1d(np.append(waves_native, waves_native[-1]*2 - waves_native[-2]),
                   np.append(P, P[-1]), kind = 'nearest')

    model = modelfn(waves = waves, fluxes = ifn(waves), sub_xs = sub_xs, PSF_fns = PSF_fns, dxs = dxs, dys = dys, modX = modX, modY = modY)
    fluxes = all_data[0]["fluxes"]
    dfluxes = all_data[0]["dfluxes"]

    print("*")
    pulls = ((fluxes - model)/dfluxes).flatten()
    return pulls
    
def get_per_pix(verbose = False):
    AB_25 = (0.11/waves**2.)*10.**(-0.4*25)

    fluxes2D = modelfn(waves, fluxes = AB_25*flamb_to_photons_per_wave, sub_xs = sub_xs, PSF_fns = PSF_fns,
                     dxs = dxs, dys = dys,
                     modX = modX, modY = modY)

    dfluxes2D = np.sqrt(fluxes2D + 1.25*exp_time)

    if verbose:
        save_img(fluxes2D, "data_per_pix.fits")

    P, NA, Cmat = miniLM_new(ministart = np.ones(len(waves_native), dtype=np.float64), miniscale = np.ones(len(waves_native), dtype=np.float64),
                             residfn = pullfn_per_pix, passdata = dict(fluxes = fluxes2D, dfluxes = dfluxes2D), verbose = True, maxiter = 100, save_jacobian = True)
    
    if verbose:
        save_img(Cmat, "Cmat_per_pix.fits")

    oneD_uncs = np.sqrt(np.diag(Cmat))
    print("oneD_uncs", oneD_uncs)
    StoN = P/oneD_uncs
    print("total StoN ", np.sqrt(np.dot(StoN, StoN)))
    
    return P, Cmat

def get_PSF(rel_x, rel_y, wave, linear_coeff_fns, PSF_fns):
    tot_weight = 0.
    tmp_mod = 0.
    
    for wave_key in PSF_fns:
        weight = linear_coeff_fns[wave_key](wave)
        tot_weight += weight
        
        if weight > 0:
            tmp_mod += PSF_fns[wave_key](rel_y, rel_x)*weight#rel_x, rel_y)*weight
            #print("wave", wave, wave_key, weight)
    assert np.isclose(tot_weight, 1.)
    return tmp_mod

class prism:
    def __init__(self, total_exp_time, dxs = [0.2, 0.7, 0.45, 0.95], dys = [0.2, 0.2, 0.7, 0.7], LSF_scale = 1.0, R_scale = 1.0):
        self.subsample = 5.
        self.total_exp_time = total_exp_time # Over all dithers
        
        self.dxs = dxs
        self.dys = dys
        """
        for dx in np.arange(2, dtype=np.float64)/2.:
            for dy in np.arange(2, dtype=np.float64)/2.:
                self.dxs.append(dx + 0.2)
                self.dys.append(dy + 0.2)
        """

        self.exp_time = self.total_exp_time/len(self.dxs)

        self.waves_native, self.dwaves_native = resolution_to_wavelengths(source_dir = wfirst_data_path + "/pixel-level/input/", IFURfl = "Prism_R.txt",
                                                                          min_wave = 7500., max_wave = 18000., R_scale = R_scale)

        self.wave_to_x_fn = interp1d(self.waves_native, np.arange(len(self.waves_native), dtype=np.float64), kind = 'linear')

        self.waves = interp1d(np.arange(len(self.waves_native), dtype=np.float64), self.waves_native, kind = 'linear')(
            np.arange(0, len(self.waves_native) - 1., 1./self.subsample, dtype=np.float64))

        self.dwaves = self.waves[1:] - self.waves[:-1]
        self.dwaves = np.append(self.dwaves, self.dwaves[-1])

        self.sub_xs = self.wave_to_x_fn(self.waves)

        self.eff_area = file_to_fn(wfirst_data_path + "/pixel-level/input/P100.txt", kind = 'linear')
        self.flamb_to_photons_per_wave = flamb_to_photons_per_wave(flamb = np.ones(len(self.waves), dtype=np.float64),
                                                                   meters2 = self.eff_area(self.waves), waves = self.waves, dwaves = self.dwaves)*self.exp_time
        self.photons_per_wave_to_flamb = photons_per_wave_to_flamb(photons_per_wave = np.ones(len(self.waves), dtype=np.float64)/self.exp_time,
                                                                   meters2 = self.eff_area(self.waves), waves = self.waves, dwaves = self.dwaves)

        self.PSFs = initialize_PSFs(pixel_scales = [22], slice_scales = [22],
                                    PSF_source = "Prism_3Element_PSF_SCA08",
                                    path = wfirst_data_path + "/pixel-level/") # Has pixel convolution

        self.psf_x = np.arange(len(self.PSFs[(8000, 22, 22)]), dtype=np.float64)/22.
        self.psf_x -= np.mean(self.psf_x)

        self.PSF_fns = {}
        self.linear_coeff_fns = {}

        all_PSFs = []

        for i, wave in enumerate(self.PSFs["waves"]):
            self.PSFs[(wave, 22, 22)][0,:] = 0.
            self.PSFs[(wave, 22, 22)][-1,:] = 0.
            self.PSFs[(wave, 22, 22)][:,0] = 0.
            self.PSFs[(wave, 22, 22)][:,-1] = 0.

            self.PSF_fns[wave] = RectBivariateSpline(self.psf_x, self.psf_x*LSF_scale, self.PSFs[(wave, 22, 22)]/LSF_scale,
                                                     kx = 1, ky = 1) #interp2d(psf_x, psf_x, PSFs[(wave, 22, 22)], fill_value = 0.)
            y_vals = np.zeros(len(self.PSFs["waves"]))
            y_vals[i] = 1.

            self.linear_coeff_fns[wave] = interp1d(self.PSFs["waves"], y_vals, kind = 'linear')
            all_PSFs.append(self.PSFs[(wave, 22, 22)])
            all_PSFs[-1][0,0] = wave

        #save_img(all_PSFs, "all_PSFs.fits")
        #modX, modY = np.meshgrid(np.arange(len(waves) + 10, dtype=np.float64), np.arange(15, dtype=np.float64))
        self.modX = np.arange(len(self.waves_native) + 10, dtype=np.float64)
        self.modX -= 5

        self.modY = np.arange(15, dtype=np.float64)
        self.modY -= np.mean(self.modY)
    

    def modelfn(self, flamb):
        fluxes = flamb*self.flamb_to_photons_per_wave

        model = []
        for j in range(len(self.dxs)):
            tmp_mod = np.zeros([len(self.modY), len(self.modX)], dtype=np.float64)

            for i in range(len(self.waves)):
                tmp_mod += fluxes[i]*get_PSF(self.modX - (self.sub_xs[i] + self.dxs[j]), self.modY - self.dys[j], wave = self.waves[i], linear_coeff_fns = self.linear_coeff_fns, PSF_fns = self.PSF_fns)
            model.append(tmp_mod)
        model = np.array(model)

        return model


if __name__ == "__main__":
    the_prism = prism()
    #get_per_pix(the_prism, verbose = True)
