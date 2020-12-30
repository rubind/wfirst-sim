import numpy as np
import sncosmo
import os
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
from tqdm import trange


def double_gauss(mu, sigp, sigm, size):
    """Double Gaussian distribution. Note: mu is the mode and not the mean."""
    
    sigm = abs(sigm) # Just in case

    p = np.array([sigp, sigm], dtype=np.float64) # probability of one side is proportional to sigma on that side
    p /= sum(p)


    sig = np.random.choice([sigp, -sigm], size = size, replace = True, p = p)
    
    return abs(np.random.normal(size = size))*sig + mu



def make_SALT2_params(size):
    """Generates "noise"-free mB, x1, c. Trained on JLA SALT2-4 SDSS z < 0.2 and SNLS z < 0.5. Very simple model with linear alpha/beta and same distribution irrspective of host-mass. mB needs h=0.7 distance modulus added to it."""

    color = double_gauss(-0.0474801042369, 0.0965032273527, 0.0428443663595, size = size)
    x1 = double_gauss(0.872727291354, 0.358731835038, 1.42806797468, size = size)
    mass = double_gauss(10.701690617, 0.334359086569, 1.0750402101, size = size)

    mB = -19.0199168813 - 0.0838387899933*(mass > 10.) + 3.20907949118*color - 0.137042055737*x1

    return mB, x1, color, mass






def get_SNCosmo_model(redshift, x1, c, MV, daymax, source):
    #sncosmo_model = sncosmo.Model(source="salt2-extended")
    sncosmo_model = sncosmo.Model(source=source)

    cosmo = FlatLambdaCDM(Om0 = 0.3, H0 = 70.)
    ten_pc_z = 2.33494867e-9
    assert abs(cosmo.distmod(z=ten_pc_z).value) < 1.e-3, "Distance modulus zeropoint wrong!"

    ampl = 1.
    for i in range(2):
        sncosmo_model.set(z=ten_pc_z, t0=0., x0=ampl, x1 = x1, c = c)

        mag = sncosmo_model.bandmag('bessellv', 'ab', 0.)
        #print "mag, ampl ", mag, ampl
        ampl *= 10.**(0.4*(mag - MV))
        
    mu = cosmo.distmod(redshift).value
    #print "mu ", mu

    sncosmo_model.set(z=redshift, t0=daymax, x0 = ampl*10**(-0.4*mu))
    return sncosmo_model



def file_to_fn(fl, col = 1):
    vals = np.loadtxt(fl)
    x = vals[:,0]
    y = vals[:,col]

    return interp1d(x, y, kind = 'linear', bounds_error = False, fill_value = 0.)

def get_SNe_per_bin(square_degrees, survey_duration, max_redshift, redshift_step):
    wfirst_data_path = os.environ["WFIRST_SIM_DATA"]
    rates_fn_full = file_to_fn(wfirst_data_path + "/pixel-level/input/SN_rates.txt") # per 1e-4 year Mpc^3 (h = 0.7)
    
    redshift_set = np.arange(redshift_step/2., max_redshift + redshift_step/2., redshift_step)
    print("redshift_set", redshift_set)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    frac_of_sky = (square_degrees)/(4.*np.pi*(180./np.pi)**2.)
    print("frac_of_sky", frac_of_sky)

    volume_in_shells = cosmo.comoving_volume(redshift_step/2. + redshift_set).value - cosmo.comoving_volume(redshift_set - redshift_step/2.).value
    SNe_in_survey_field = volume_in_shells*frac_of_sky*rates_fn_full(redshift_set) * survey_duration/(1. + redshift_set) * 1e-4

    print (SNe_in_survey_field, sum(SNe_in_survey_field))
    for item in zip(redshift_set, SNe_in_survey_field):
        print(item)
    return redshift_set, SNe_in_survey_field


def make_redshifts(square_degrees, survey_duration, max_redshift, redshift_step = 0.05):
    redshift_set, SNe_in_survey_field = get_SNe_per_bin(square_degrees = square_degrees, survey_duration = survey_duration, max_redshift = max_redshift, redshift_step = redshift_step)

    SNe_actual = [np.random.poisson(item) for item in SNe_in_survey_field]
    redshifts = []
    for i in range(len(redshift_set)):
        for j in range(SNe_actual[i]):
            redshifts.append(redshift_set[i])
    return np.array(redshifts)


def realize_SNe(redshifts, daymaxes, source):
    MBs, x1s, colors, masses = make_SALT2_params(size = len(redshifts))
    HRs = np.random.normal(size = len(redshifts))*0.1 + 0.055*redshifts*np.random.normal(size = len(redshifts)) + 5/np.log(10.)*(0.001/redshifts)*np.random.normal(size = len(redshifts))

    MBs += HRs
    
    sncosmo_models = []
    for i in trange(len(redshifts)):
        sncosmo_model = get_SNCosmo_model(redshift = redshifts[i],
                                          x1 = x1s[i],
                                          c = colors[i],
                                          MV = MBs[i] - colors[i],
                                          daymax = daymaxes[i], source = source)
        sncosmo_models.append(sncosmo_model)

    return MBs - colors, x1s, colors, masses, sncosmo_models, HRs


def get_SN_population(square_degrees = None, survey_duration = None, max_redshift = None, redshift_step = 0.05, scramble_z_order = 0, redshifts = None):
    if square_degrees != None:
        redshifts = make_redshifts(square_degrees = square_degrees, survey_duration = survey_duration, max_redshift = max_redshift, redshift_step = redshift_step)

        if scramble_z_order:
            np.random.shuffle(redshifts)

    wfirst_data_path = os.environ["WFIRST_SIM_DATA"]
    source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")
    
    daymaxes = np.random.random(size = len(redshifts))*survey_duration*365.24
    MVs, x1s, colors, masses, sncosmo_models, HRs = realize_SNe(redshifts, daymaxes, source)
    return redshifts, daymaxes, MVs, x1s, colors, masses, sncosmo_models, HRs


# import os
# import sys
# wfirst_path = os.environ["WFIRST"]
# sys.path.append(wfirst_path + "/scripts/stan_cosmo/")
# source = sncosmo.SALT2Source(modeldir=wfirst_data_path + "/salt2_extended/")


# get_SN_population(square_degrees = 1., survey_duration = 2., max_redshift = 2.5, redshift_step = 0.05)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    mB, x1, color, mass = make_SALT2_params(100000)

    print ("Median mB", np.median(mB))


    plt.subplot(2,2,1)
    plt.hist(mB, bins = 30)
    plt.title("mB")

    plt.subplot(2,2,2)
    plt.hist(x1, bins = 30)
    plt.title("x1")

    plt.subplot(2,2,3)
    plt.hist(color, bins = 30)
    plt.title("color")

    plt.subplot(2,2,4)
    plt.hist(mass, bins = 30)
    plt.title("mass")
    
    plt.savefig("SALT_params.pdf")
    plt.close()


    redshifts = make_redshifts(square_degrees = 100., survey_duration = 1., max_redshift = 2.0, redshift_step = 0.05)
    plt.hist(redshifts, bins = np.arange(0., 2.01, 0.1))
    plt.savefig("redshifts.pdf")
    plt.close()
