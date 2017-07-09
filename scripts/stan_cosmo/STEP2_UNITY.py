from numpy import *
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import pystan
import time
import cPickle as pickle
from scipy.interpolate import interp1d
import pyfits
import sys
import os
wfirst_path = os.environ["WFIRST"]
sys.path.append(wfirst_path + "/scripts/pixel-level/")
from pixel_level_ETC2 import initialize_PSFs, get_spec_with_err, get_sncosmo, save_img, spectrum_to_matched_resolution, bin_vals_fixed_bins, interpfile
import commands
import argparse
import extinction


wfirst_data_path = os.environ["WFIRST_SIM_DATA"]

def get_sys_coeff_fns(xmin, xmax, xscale, ncoeff, scale, xsteps = 50, epsilon = 1e-6):
    assert xsteps > ncoeff*2, "Run with more xsteps!"

    xvals = linspace(xmin - epsilon, xmax + epsilon, xsteps)

    cmat = zeros([xsteps]*2, dtype=float64)
    for i in range(xsteps):
        for j in range(i, xsteps):
            cmat[i,j] = (1. + (xvals[i] - xvals[j])**2. / (2. * xscale**2.) )**(-1.)
            cmat[j,i] = cmat[i,j]
    cmat += diag([1e-6]*xsteps)

    evals, evecs = linalg.eig(cmat)

    evals = real(evals)
    evecs = real(evecs)
    evecs = transpose(evecs)

    inds = argsort(evals)[::-1]

    evals = evals[inds]
    evecs = evecs[inds]

    for i in range(xsteps):
        evecs[i] *= sqrt(evals[i])

    if median(evecs[0]) < 0:
        print "Flipping sign!"
        evecs *= -1

    """
    save_img(cmat, "cmat.fits")
    new_mat = cmat*0.

    for i in range(xsteps):
        new_mat += outer(evecs[i], evecs[i])

    save_img(new_mat, "new_mat.fits")
    """
    
    return [interp1d(xvals, evecs[i]*scale, kind = 'linear', bounds_error = False, fill_value = NaN) for i in range(ncoeff)]




"""
coeffs = get_sys_coeff_fns(xmin = -3, xmax = 4, xscale = 1, ncoeff = 15, scale = 0.005, xsteps = 50, epsilon = 1e-6)
for coeff in coeffs:
    xvals = linspace(-3, 4, 50)
    plt.plot(xvals, coeff(xvals))
plt.savefig("CRNL_coeffs.pdf")

fsdfadsf

"""
"""
def get_sys_coeff_fns_spline(xmin, xmax, ncoeff, xsteps = 50):
    assert xsteps > ncoeff*2, "Run with more xsteps!"

    xvals = linspace(xmin, xmax, xsteps)

    cmat = zeros([xsteps]*2, dtype=float64)

    for i in range(5):
        xnodes = linspace(xmin, xmax, 5)
        ynodes = xnodes*0
        ynodes[i] = 1

        ifn = interp1d(xnodes, ynodes, kind = 'linear')
        cmat += outer(ifn(xvals), ifn(xvals))

    cmat += diag([1e-6]*xsteps)

    evals, evecs = linalg.eig(cmat)

    evals = real(evals)
    evecs = real(evecs)
    evecs = transpose(evecs)

    inds = argsort(evals)[::-1]

    evals = evals[inds]
    evecs = evecs[inds]

    for i in range(xsteps):
        evecs[i] *= sqrt(evals[i])

    save_img(cmat, "cmat.fits")
    new_mat = cmat*0.

    for i in range(xsteps):
        new_mat += outer(evecs[i], evecs[i])

    save_img(new_mat, "new_mat.fits")
    
    for i in range(ncoeff):
        plt.plot(xvals, evecs[i])
    #print evals[:100]
    #print [evals[:i+1].sum()/evals.sum() for i in range(0,300, 50)]

"""



def get_eigen():
    # Results from a very simple at-max SNf PCA. This is something we're still working on!
    fl = pyfits.open(wfirst_data_path + "/pixel-level/input/eigen_vects_norm.fits")
    evecs = fl[0].data
    fl.close()
    
    eigen_fns = []
    for i in range(1,1+opts.neigen):
        eigen_fns.append(interp1d(evecs[0], evecs[i], kind = 'linear', fill_value = 0, bounds_error = False))
    return eigen_fns


def getIGextinction(z, efflamb_A):
    efflamb = efflamb_A/10000.

    restA = z*(0.12835272772310144 - 0.5549447330170015*efflamb + 1.2280496123038864*efflamb**2 - 1.5826072252730259*efflamb**3 + 1.2123528947173046*efflamb**4 - 0.5400476644358357*efflamb**5 + 0.12868904067282008*efflamb**6 - 0.012655137320257642*efflamb**7 + z*(-0.017751821448063823 + 0.07672233425134578*efflamb - 0.18766993140145194*efflamb**2 + 0.2705543905364898*efflamb**3 - 0.22705950051851762*efflamb**4 + 0.1082170905126371*efflamb**5 - 0.02708449897611635*efflamb**6 + 0.0027608766952176883*efflamb**7))
    
    return restA


def get_jacobian(rest_waves, redshift, params, electrons_per_sec, fund_cal_fns, crnl_fns):

    # "PSF_wghtd_e_allbutdark"

    if redshift > 0.1:
        EBV_step = 0.01*params["mw_norm"]
    else:
        EBV_step = 0.05*params["mw_norm"]

    jacobian = []
    jacobian_names = []

    for fund_cal_fn in fund_cal_fns:
        jacobian.append(fund_cal_fn(rest_waves*(1. + redshift)))
        jacobian_names.append("Fundamental Calibration")

    for crnl_fn in crnl_fns:
        jacobian.append(crnl_fn(electrons_per_sec))
        jacobian_names.append("Count-Rate Nonlinearity")

    ff = extinction.Fitzpatrick99(3.1)
    jacobian.append(ff(rest_waves*(1. + redshift), 3.1*EBV_step))
    jacobian_names.append("MW Extinction Normalization")

    jacobian.append(ff(rest_waves*(1. + redshift), 3.1*params["mw_zeropoint"])) # Independent zeropoint uncertainty of this in E(B-V)
    jacobian_names.append("MW Extinction Zeropoint")

    if redshift < 0.1:  # R_V uncertainty of 0.2; E(B-V) of 0.01; ignore nearby SNe as these are spread on sky.
        jacobian.append([0]*rest_waves)
    else:
        R_V_prime = 3.1 - params["mw_RV_uncertainty"]
        ff_prime = extinction.Fitzpatrick99(R_V_prime)

        jacobian.append(ff_prime(rest_waves*(1. + redshift), 3.1*EBV_step) - ff(rest_waves*(1. + redshift), 3.1*EBV_step))
    jacobian_names.append("MW Extinction $R_V$")

    
    jacobian.append(getIGextinction(redshift, rest_waves)*params["ig_extinction"])
    jacobian_names.append("IG Extinction")
    

    jacobian = transpose(array(jacobian))
    #save_img(jacobian, "jacobian.fits")
    
    return jacobian, jacobian_names


def get_F99(this_rest_lambs):
    ff = extinction.Fitzpatrick99(3.1)
    fprime = extinction.Fitzpatrick99(3.11)
    F99_31 = ff(this_rest_lambs, 3.1)
    F99_dAdRV = (fprime(this_rest_lambs, 3.11) - F99_31)/0.01
    
    return F99_31, F99_dAdRV


def generate_data(SN_data):
    print "Reading PSFs..."

    PSFs = initialize_PSFs(pixel_scales = [10, 15, 30], slice_scales = [30, 30, 30], PSF_source = SN_data["survey_parameters"]["PSFs"]) # Native in units of 5 mas, so this is 0".075 and 0".15

    args = {"pixel_scale": SN_data["survey_parameters"]["IFU_pixel_scale"],
            "slice_scale": SN_data["survey_parameters"]["IFU_slice_in_pixels"]*SN_data["survey_parameters"]["IFU_pixel_scale"],
            "offset_par": int(around(SN_data["survey_parameters"]["IFU_pixel_scale"]*0.25/0.005)), "offset_perp": 0,
            "dark_current": opts.IFCdark, "mdl": 'hsiao',
            "PSFs": PSFs, "IFURfl": SN_data["survey_parameters"]["IFU_resolution"],
            "IPC": opts.IFCIPC, "bad_pixel_rate": SN_data["survey_parameters"]["bad_pixel_rate"],
            "min_wave": opts.IFCminwave,
            "max_wave": opts.IFCmaxwave,
            "read_noise_floor": opts.IFCRNfloor,
            "effareafl": SN_data["survey_parameters"]["IFU_effective_area"],
            "TTel": opts.TTel,
            "source_dir": wfirst_data_path + "/pixel-level/input/"}

    

    ETC_result = get_spec_with_err(redshift = 0.1, exp_time = 1, show_plots = 0, phase = 0, **args)
    args["waves"] = ETC_result["obs_waves"]
    high_res_obs_waves = arange(ETC_result["obs_waves"].min() - 10., ETC_result["obs_waves"].max() + 11., 10.)

    if not all(args["waves"] == SN_data["IFC_waves"]):
        
        for j in range(len(SN_data["SN_observations"])):
            SN_data["SN_observations"][j]["gal_background"] = interp1d(SN_data["IFC_waves"], SN_data["SN_observations"][j]["gal_background"], kind = 'linear', bounds_error = False, fill_value = median(SN_data["SN_observations"][j]["gal_background"]))(args["waves"])
        SN_data["IFC_waves"] = args["waves"]
        

    args["PSFs"] = ETC_result["PSFs"]
    

    print "Setting up..."

    eigen_fns = get_eigen()

    has_IFS_mask = zeros(len(SN_data["SN_observations"]))

    for i in range(len(SN_data["SN_observations"])):
        if len(SN_data["SN_observations"][i]["IFS_dates"]) > 0:
            phases = [(IFS_date - SN_data["SN_table"]["daymaxes"][i])/(1. + SN_data["SN_table"]["redshifts"][i]) for IFS_date in SN_data["SN_observations"][i]["IFS_dates"]]
            phases = array(phases)
            if any(abs(phases) < 5):
                has_IFS_mask[i] = 1
            else:
                print "IFS observations found, but not with 5 rest-frame days of max!", phases


    has_IFS_inds = where(has_IFS_mask)[0]
    inds_in_order_of_z = argsort(array(SN_data["SN_table"]["redshifts"][has_IFS_inds]))
    has_IFS_inds = has_IFS_inds[inds_in_order_of_z]

    assert sum(has_IFS_inds) == sum(where(has_IFS_mask)[0])

    print has_IFS_mask
    print len(has_IFS_mask)
    
    

    redshifts = sort(unique(array(SN_data["SN_table"]["redshifts"][has_IFS_inds])))
    redshifts = concatenate(([0.05], redshifts))
    print "redshifts ", redshifts

    
    nred = len(redshifts)

    print "neigen", opts.neigen
    
    nsyscoeff = 10

    nsys = (nsyscoeff + # Fundamental
            nsyscoeff + # CRNL
            4) # MW norm, zeropoint, RV, IG extinction

    
    params = dict(mw_norm = opts.mwnorm, mw_zeropoint = opts.mwZP, mw_RV_uncertainty = opts.mwRV, ig_extinction = opts.IGext, crnl = opts.crnl, fund_cal = opts.fund)

    crnl_fns = get_sys_coeff_fns(xmin = -3, xmax = 1., xscale = 1., ncoeff = nsyscoeff, scale = params["crnl"], xsteps = 50)
    fund_cal_fns = get_sys_coeff_fns(xmin = 2999.9, xmax = 21000.1, xscale = 5000., ncoeff = nsyscoeff, scale = params["fund_cal"], xsteps = 50)
    if opts.nredcoeff > 2:
        scale_factor_fns = get_sys_coeff_fns(xmin = 1./(1. + 2.0), xmax = 1., xscale = 0.5, ncoeff = opts.nredcoeff, scale = 1.0, xsteps = 50)
    elif opts.nredcoeff == 2:
        scale_factor_fns = [lambda x: 1., lambda x: (x - 1.)]
    elif opts.nredcoeff == 1:
        scale_factor_fns = [lambda x: 1.]
    else:
        raise Exception("Other number of nredcoeff! " + str(opts.nredcoeff))

    rest_waves = exp(linspace(log(opts.bluewave), log(opts.redwave), opts.nrestlamb))
    
    spectral_resolution = (log(opts.redwave) - log(opts.bluewave))/(opts.nrestlamb/2.)

    rest_waves_pad = concatenate(([rest_waves[0]/(rest_waves[1]/rest_waves[0])], rest_waves, [rest_waves[-1]*(rest_waves[-1]/rest_waves[-2])]))
    print "rest_waves_pad", rest_waves_pad


    scale_factors = 1./(1 + redshifts)

    redshift_coeffs = zeros([nred, opts.nredcoeff], dtype=float64)

    for i in range(opts.nredcoeff):
        redshift_coeffs[:, i] = scale_factor_fns[i](scale_factors)
    plt.plot(redshifts, redshift_coeffs)
    plt.savefig("redshift_coeffs.pdf")
    plt.close()
    
    zinds = [0]*opts.nnearby
    for ind in has_IFS_inds:
        zinds.append(list(redshifts).index(SN_data["SN_table"]["redshifts"][ind]))

    plt.plot(zinds)
    plt.savefig("zinds.pdf")
    plt.close()

    nsne = len(zinds)
    print "nsne ", nsne

    true_EBVs = random.exponential(size = nsne)*0.1
    redshift_vector = array([redshifts[zind] for zind in zinds])

    true_mags = random.normal(size = nsne)*sqrt(opts.gray**2. + 0.055**2. * redshift_vector**2.  + 0.00217147**2. / redshift_vector**2. )
    true_RVs = random.normal(size = nsne)*0.31 + 3.1


    rest_mod, NA, NA, NA = get_sncosmo('hsiao', 1., rest_waves*2., array([1e4]), phase = 0.)
    rest_mod /= sum(rest_mod) # This is just for intializing the model; it's not used anywhere else.

    true_projs = random.normal(size = (nsne, opts.neigen))


    true_fluxes = []
    fluxes = []
    dfluxes = []

    dflux_dsys = zeros((nsne, opts.nrestlamb, nsys), dtype=float64)


    for i in range(nsne):
        this_redshift = redshifts[zinds[i]]

        print "*"*20 + " this_redshift ", this_redshift

        

        if this_redshift > 0.1:
            this_rest_lambs = args["waves"]/(1. + this_redshift)
        else:
            this_rest_lambs = exp(arange(log(3200.), log(9000.), 0.005))


        if this_redshift > 0.1:
            # WFIRST
            f_lamb, NA, NA, NA = get_sncosmo('hsiao', this_redshift, high_res_obs_waves, array([1e4]), phase = 0., absmag = -19.08 - 3.1*median(true_EBVs))
            this_eigen = array([item(high_res_obs_waves/(1 + this_redshift)) for item in eigen_fns])
            F99_31, F99_dAdRV = get_F99(high_res_obs_waves/(1 + this_redshift))

        else:
            # Nearby
            f_lamb, NA, NA, NA = get_sncosmo('hsiao', this_redshift, this_rest_lambs*(1 + this_redshift), array([1e4]), phase = 0., absmag = -19.08 - 3.1*median(true_EBVs))
            this_eigen = array([item(this_rest_lambs) for item in eigen_fns])
            F99_31, F99_dAdRV = get_F99(this_rest_lambs)
        


        f_lamb *= 10**(-0.4*(true_mags[i] + true_EBVs[i]*((true_RVs[i] - 3.1)*F99_dAdRV + F99_31) + dot(true_projs[i], this_eigen)
                             ))


        IFUR_FN = interpfile(args["source_dir"] + "/" + args["IFURfl"])

        if this_redshift > 0.1:
            at_max_mdl = spectrum_to_matched_resolution(high_res_obs_waves, f_lamb, min_wave = args["min_wave"], max_wave = args["max_wave"], IFUR = IFUR_FN, pixel_scale = args["pixel_scale"])
            args["mdl"] = at_max_mdl

            orig_ind = has_IFS_inds[i - opts.nnearby]

            assert len(SN_data["SN_observations"][orig_ind]["IFS_dates"]) > 0
 
            total_spec_StoN = 0.
            
            for j in range(len(SN_data["SN_observations"][orig_ind]["IFS_dates"])):
                phase = (SN_data["SN_observations"][orig_ind]["IFS_dates"][j] - SN_data["SN_table"]["daymaxes"][orig_ind])/(1. + SN_data["SN_table"]["redshifts"][orig_ind])
                if abs(phase) < 5:
                    print SN_data["SN_observations"][orig_ind]["IFS_dates"][j], SN_data["SN_table"]["daymaxes"][orig_ind], SN_data["SN_table"]["redshifts"][orig_ind], phase
                
                    this_ETC_result = get_spec_with_err(redshift = 0., exp_time = SN_data["SN_observations"][orig_ind]["IFS_exptimes"][j], show_plots = 0,
                                                        gal_flamb = interp1d(SN_data["IFC_waves"], SN_data["SN_observations"][j]["gal_background"], kind = 'linear'), phase = 0, **args)
                    total_spec_StoN += this_ETC_result["spec_S/N"]**2.
                    print "S/N", j, median(this_ETC_result["spec_S/N"])
                    exp_time = SN_data["SN_observations"][orig_ind]["IFS_exptimes"][j]

            total_spec_StoN = sqrt(total_spec_StoN)
            print "Total S/N, no ref", median(total_spec_StoN)
            
            ETC_result = {"f_lamb_SN": at_max_mdl, "PSF_wghtd_e_allbutdark": this_ETC_result["PSF_wghtd_e_allbutdark"], "exp_time": exp_time}
            noise_in_SN_active = ETC_result["f_lamb_SN"]/total_spec_StoN

            args["mdl"] = at_max_mdl/100. # For the reference, the SN is faint
            phase_ref = (SN_data["SN_observations"][orig_ind]["IFS_dates"][-1] - SN_data["SN_table"]["daymaxes"][orig_ind])/(1. + SN_data["SN_table"]["redshifts"][orig_ind])
            assert phase_ref > 60, "No reference, or it's not the last observation!"
            


            ref_SN = 0
            total_refs_found = 0
            

            for j in range(len(SN_data["SN_observations"][orig_ind]["IFS_dates"])):
                phase = (SN_data["SN_observations"][orig_ind]["IFS_dates"][j] - SN_data["SN_table"]["daymaxes"][orig_ind])/(1. + SN_data["SN_table"]["redshifts"][orig_ind])
                if phase > 60:
                    ETC_result_ref = get_spec_with_err(redshift = 0., exp_time = SN_data["SN_observations"][orig_ind]["IFS_exptimes"][j], show_plots = 0, phase = 0,
                                                       gal_flamb = interp1d(SN_data["IFC_waves"], SN_data["SN_observations"][j]["gal_background"], kind = 'linear'), **args)
                    ref_SN += ETC_result_ref["spec_S/N"]**2.
                    total_refs_found += 1

            assert total_refs_found == SN_data["survey_parameters"]["number_of_reference_dithers"], "Found %i refs, but number_of_reference_dithers %i %s" % (total_refs_found, SN_data["survey_parameters"]["number_of_reference_dithers"], str(SN_data["SN_observations"][orig_ind]["IFS_dates"]))

            ref_SN = sqrt(ref_SN)
            noise_in_reference = args["mdl"]/ref_SN

            total_noise = sqrt(noise_in_SN_active**2 + noise_in_reference**2)
            total_spec_StoN = at_max_mdl/total_noise
            
            ETC_result["spec_S/N"] = total_spec_StoN
            print "Total S/N, with ref", median(total_spec_StoN)

            
            
        else:
            resl = median(0.5*(this_rest_lambs[1:] + this_rest_lambs[:-1])/(2*(this_rest_lambs[1:] - this_rest_lambs[:-1])))

            print resl
            print len(f_lamb), len(this_rest_lambs)

            ETC_result = {"f_lamb_SN": f_lamb, "PSF_wghtd_e_allbutdark": [NaN]*len(f_lamb), "spec_S/N": ones(len(f_lamb), dtype=float64)*sqrt(resl/150.)/0.12, "exp_time": 1.}


        print "Binning ", len(this_rest_lambs), len(ETC_result["f_lamb_SN"]), len(ETC_result["spec_S/N"]), len(rest_waves_pad)
        
        bin_flux, bin_err = bin_vals_fixed_bins(this_rest_lambs, ETC_result["f_lamb_SN"] + random.normal(size = len(this_rest_lambs))*ETC_result["f_lamb_SN"]/ETC_result["spec_S/N"],
                                                ETC_result["f_lamb_SN"]/ETC_result["spec_S/N"], rest_waves_pad)
        bin_flux = bin_flux[1:-1] # Eliminate the padding
        bin_err = bin_err[1:-1] # Eliminate the padding
        
        
        bin_e_per_sec, NA = bin_vals_fixed_bins(this_rest_lambs, log10(array(ETC_result["PSF_wghtd_e_allbutdark"])/ETC_result["exp_time"]), ones(len(this_rest_lambs), dtype=float64), rest_waves_pad)
        bin_e_per_sec = bin_e_per_sec[1:-1]



        jacobian, jacobian_names = get_jacobian(rest_waves, redshifts[zinds[i]], params, bin_e_per_sec, fund_cal_fns, crnl_fns)
        
        inds = where(isnan(jacobian))
        jacobian[inds] = 0
        
        
        inds = where(isnan(bin_flux))
        bin_flux[inds] = 0
        bin_err[inds] = max(bin_flux)*100
            

        dflux_dsys[i] = jacobian
        fluxes.append(bin_flux)
        dfluxes.append(bin_err)


    fluxes = array(fluxes)*1e16
    dfluxes = array(dfluxes)*1e16

    save_img(fluxes, "fluxes.fits")
    save_img(dfluxes, "dfluxes.fits")
    save_img(fluxes/dfluxes, "SNRs.fits")
    save_img(dflux_dsys, "dflux_dsys.fits")

    F99_31, F99_dAdRV = get_F99(rest_waves)

    plt.plot(rest_waves, F99_31)
    plt.plot(rest_waves, F99_dAdRV)
    plt.savefig("F99.pdf")
    plt.close()



    for i in range(nsne):
        ratio = median(fluxes[i,opts.nrestlamb/4:opts.nrestlamb/2])/median(fluxes[i,opts.nrestlamb/2:3*opts.nrestlamb/4])
        plt.plot(true_EBVs[i], ratio, '.', color = 'k')
    plt.savefig("EBV_check.pdf")
    plt.close()



    eigen_vecs = array([item(rest_waves) for item in eigen_fns])

    save_img(eigen_vecs, "sampled_eigen_vecs.fits")

    stan_data = dict(nrestlamb = opts.nrestlamb, nsne = nsne, nred = nred, neigen = opts.neigen, nsys = nsys, ncoeff = opts.nredcoeff,
                     fluxes = fluxes, fluxerrs = dfluxes, eigen_vecs = eigen_vecs, dflux_dsys = dflux_dsys, CCM_31 = F99_31, CCM_dAdRV = F99_dAdRV,
                     redshift_coeffs = redshift_coeffs, zinds = zinds,
                     gray_variance = 0.055**2. * redshift_vector**2.  + 0.00217147**2. / redshift_vector**2.
                )

    other_inputs = dict(true_fluxes = true_fluxes, true_projs = true_projs, true_mags = true_mags, true_EBVs = true_EBVs, true_RVs = true_RVs, rest_mod = rest_mod, jacobian_names = jacobian_names, redshifts = redshifts)
    return stan_data, other_inputs


def initfn():
    nsne = stan_data["nsne"]
    nred = stan_data["nred"]
    neigen = stan_data["neigen"]
    nsys = stan_data["nsys"]
    ncoeff = stan_data["ncoeff"]

    rest_mod_init = other_inputs["rest_mod"]*exp(random.normal(size = stan_data["nrestlamb"])*0.1)
    rest_mod_init /= sum(rest_mod_init)

    temp_fluxes = array([stan_data["fluxes"][i]*10**(0.4*stan_data["CCM_31"]*other_inputs["true_EBVs"][i]) for i in range(nsne)])

    init_mus = []
    for i in range(0,nred):
        inds = where(array(stan_data["zinds"]) == i)
        assert len(inds[0]) > 0, "Couldn't find " + str(i) # Just a check, should never fail
        init_mus.append(-2.5*log10(median(median(temp_fluxes[inds], axis = 0)/rest_mod_init)))
    print "init_mus ", init_mus
    

    return dict(
        rest_mod = rest_mod_init,

        mags = other_inputs["true_mags"] + [init_mus[ind] for ind in stan_data["zinds"]],
        mus = init_mus,
        log_gray_disp = random.normal()*0.1 - 2.3,

        #log_true_EBV = log(other_inputs["true_EBVs"]),
        true_EBV = other_inputs["true_EBVs"],
        true_projs = other_inputs["true_projs"],

        #true_dRV = other_inputs["true_RVs"] - 3.1 + random.normal(size = nsne)*0.05,
        log_R_EBV_coeff = [-2.3 + random.normal()*0.1] + [0]*(ncoeff - 1),
        #log_EBV_star_coeff = [-2.3 + random.normal()*0.2] + [0]*(ncoeff - 1),
        #d_log_EBV_star_d_proj = random.normal(size = neigen)*0.05,
        RV_coeff = [3.1 + random.normal()*0.1] + [0]*(ncoeff - 1),

        log_R_proj_coeff = random.normal(size = (ncoeff, neigen))*0.1,
        proj_star_variable_coeff = random.normal(size = (ncoeff-1, neigen))*0.1,

        dsys = random.normal(size = nsys)*0.1
    )


print "Here are the inputs:"
print " ".join(sys.argv)
parser = argparse.ArgumentParser()

parser.add_argument('-pickle', help='pickle file to read', type=str)
parser.add_argument('-niter', help='number of iterations', type=int, default = 1000)
parser.add_argument('-nchains', help='number of chains', type=int, default = 8)
parser.add_argument('-writecsv', help='write csv files for chains', type=int, default = 1)
parser.add_argument('-nnearby', help='number of nearby SNe', type=int, default = 800)
parser.add_argument('-nrestlamb', help='number of rest-frame wavelengths', type=int, default=100)
parser.add_argument('-nredcoeff', help='number of redshift coeffs', type=int, default=3)
parser.add_argument('-neigen', help='number of eigen vectors', type=int, default=5)
parser.add_argument('-bluewave', help='blue rest wavelength limit', type=float, default=3300.)
parser.add_argument('-redwave', help='red rest wavelength limit', type=float, default=8600.)
parser.add_argument('-stan', help='location of stan code', type=str, default=wfirst_path + "/scripts/stan_cosmo/stan_code_SimultUNITY.txt")

parser.add_argument('-IFCdark', help='IFC dark current', type = float, default = -1)
parser.add_argument('-IFCRNfloor', help='IFC read noise floor', type = float, default = -1)
parser.add_argument('-IFCminwave', help='IFC minimum wave', type = float, default = -1)
parser.add_argument('-IFCmaxwave', help='IFC maximum wave', type = float, default = -1)
parser.add_argument('-IFCIPC', help='IFC IPC', type = float, default = -1)
parser.add_argument('-TTel', help='telescope temperature (K)', type=float, default=-1)

parser.add_argument('-mwnorm', help='MW extinction normalization uncertainty', type=float, default=0.05)
parser.add_argument('-mwZP', help='MW extinction zeropoint uncertainty', type=float, default=0.005)
parser.add_argument('-mwRV', help='MW extinction RV uncertainty', type=float, default=0.2)
parser.add_argument('-IGext', help='intergalactic extinction fractional uncertainty', type=float, default=0.25)
parser.add_argument('-crnl', help='count-rate nonlinearity (mag)', type=float, default=0.005)
parser.add_argument('-fund', help='fundamental calibration (mag)', type=float, default=0.005)
parser.add_argument('-gray', help='gray dispersion (mag)', type=float, default=0.08)


opts = parser.parse_args()

print opts






SN_data = pickle.load(open(opts.pickle, 'rb'))


assert SN_data["survey_parameters"]["adjust_each_SN_exp_time"] == 0, "Not set up for having a different exposure time for different SNe!"

for key in SN_data:
    print "SN_data ", key

for key in SN_data["survey_parameters"]:
    print "survey_parameters ", key


if opts.IFCdark == -1:
    opts.IFCdark = SN_data["survey_parameters"]["IFU_dark_current"]
if opts.IFCRNfloor == -1:
    opts.IFCRNfloor = SN_data["survey_parameters"]["IFU_read_noise_floor"]
if opts.IFCminwave == -1:
    opts.IFCminwave = SN_data["survey_parameters"]["IFU_min_wave"]
if opts.IFCmaxwave == -1:
    opts.IFCmaxwave = SN_data["survey_parameters"]["IFU_max_wave"]
if opts.IFCIPC == -1:
    opts.IFCIPC = SN_data["survey_parameters"]["interpixel_capacitance"]
if opts.TTel == -1:
    opts.TTel = SN_data["survey_parameters"]["telescope_temperature"]

print opts


stan_data, other_inputs = generate_data(SN_data)

print "Ready to start", time.asctime()

if opts.writecsv:
    sample_file = commands.getoutput("pwd") + "/samples.txt"
else:
    sample_file = None

fit = pystan.stan(file=opts.stan, data=stan_data,
                  iter=opts.niter, chains=opts.nchains, n_jobs = opts.nchains, refresh = 5, init = initfn,
                  sample_file = sample_file)

print "Done sampling", time.asctime()
fit_params = fit.extract(permuted = True)
print "Done extracting", time.asctime()

pickle.dump((fit_params, stan_data, other_inputs), open("fit_results.pickle", 'wb'))

print fit
print "Done printing", time.asctime()


print "Done!", time.asctime()

