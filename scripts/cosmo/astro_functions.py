from copy import deepcopy
from numpy import *
import pylab
#from FileRead import readcol
from scipy.interpolate import interp1d




def n_integrate(x_list, intfn, otherargs, pad_list = arange(21, dtype=float64)/10.):

    nadded = len(pad_list)
    x_list = append(x_list, pad_list) # make sure that these values are evaluated

    if len(x_list) % 2 == 0:
        nadded += 1
        x_list = append(x_list, 2.1) # make sure it's odd

    x_indices = argsort(x_list)
    x_sort_list = take(x_list, x_indices)

    #print x_sort_list

    #2n + 1 samples

    #(1/3)*(dx)*(f0 + 4*(f1 + f3 + ... + f2n-1) + 2*(f2 + f4 + ... + f2n-2) + f2n)


    intfn_sort = intfn(x_sort_list, otherargs) # (2n + 1)

    interxlist = 0.5*(x_sort_list[:-1] + x_sort_list[1:])
    
    intfn_sort_inter = intfn(interxlist, otherargs)  # (2n)


    the_integral = zeros(len(intfn_sort), dtype=float64)

    
    for i in range(1,len(intfn_sort)):
        the_integral[i] = the_integral[i - 1] + (intfn_sort[i - 1] + 4.*intfn_sort_inter[i - 1] + intfn_sort[i])*(x_sort_list[i] - x_sort_list[i - 1])/6.


    tmplist = deepcopy(the_integral)
    put(the_integral, x_indices, tmplist)

    the_integral = the_integral[:-nadded]

    return the_integral

"""
def H_inv(z, cparams):
    Om, w0, wa = cparams
    return 1./sqrt(Om*(1. + z)**3. + (1 - Om)*exp(-3.*wa*z/(1. + z))* (1. + z)**(3.*(1 + w0 + wa)))

def H_inv_Om(z, cparams):
    Om, w0, wa = cparams    
    return -((1. + z)**3 - (1 + z)**(3*(1 + w0 + wa))/exp((3*wa*z)/(1 + z)))/(2.*(Om*(1 + z)**3 + ((1 - Om)*(1 + z)**(3*(1 + w0 + wa)))/exp((3*wa*z)/(1 + z)))**1.5)

def H_inv_w0(z, cparams):
    Om, w0, wa = cparams
    return (-3*(1 - Om)*(1 + z)**(3*(1 + w0 + wa))*log(1 + z))/(2.*exp((3*wa*z)/(1 + z))*(Om*(1 + z)**3 + ((1 - Om)*(1 + z)**(3*(1 + w0 + wa)))/exp((3*wa*z)/(1 + z)))**1.5)

def H_inv_wa(z, cparams):
    Om, w0, wa = cparams
    return (3*(-1 + Om)*(1 + z)**(-1 + 3*w0 + 3*wa)*(-z + (1 + z)*log(1 + z)))/(2.*(exp((3*wa*z)/(1 + z))*Om - (-1 + Om)*(1 + z)**(3*(w0 + wa)))*sqrt(Om*(1 + z)**3 + ((1 - Om)*(1 + z)**(3*(1 + w0 + wa)))/exp((3*wa*z)/(1 + z))))
"""

def H_inv(z, cparams):
    Om, wp, wa, zp = cparams
    return 1./sqrt(Om*(1. + z)**3. + (1 - Om)*exp(-3.*wa*z/(1. + z)) * (1. + z)**(3.*(1 + wp + wa + zp + wp*zp)/(1. + zp)))

def H_inv_Om(z, cparams):
    Om, wp, wa, zp = cparams    
    return -0.5*H_inv(z, cparams)**3. * (
        (1. + z)**3. - exp(-3.*wa*z/(1. + z)) * (1. + z)**(3.*(1 + wp + wa + zp + wp*zp)/(1. + zp))
        )

def H_inv_wp(z, cparams):
    Om, wp, wa, zp = cparams
    return -0.5*H_inv(z, cparams)**3. * (
        -3. * exp(-3. * wa * z/(1. + z)) * (Om - 1.) * (1. + z)**(3.*(1 + wa + wp + zp + wp*zp)/(1. + zp)) * log(1. + z)
        )

def H_inv_wa(z, cparams):
    Om, wp, wa, zp = cparams
    return -0.5*H_inv(z, cparams)**3. * (
        (3*exp(-3*wa*z/(1. + z))*(Om - 1.)*(1 + z)**(2. + 3*wp + 3*wa/(1. + zp)) * (z*(1. + zp) - (1. + z)*log(1. + z))/
        (1 + zp))
        )

def shift_parameter_Wmat(fractional_error, zp):
    el12 = -0.0288451
    el11 = 0.105128
    el22 = 0.00791457
    el33 = 0.0005823800298267789 - 0.003937215276197356*zp + 0.008844122979160287*zp**2 - 0.007659268280406035*zp**3 + 0.0025819308752891724*zp**4
    el23 = 0.002163214336871038 - 0.007784496510622427*zp + 0.006842603532185611*zp**2 - 0.004234458343122008*zp**3 + 0.0012252910896924294*zp**4
    el13 = -0.007883962843451093 + 0.028371058844540536*zp - 0.024938274067769745*zp**2 + 0.015432734366768336*zp**3 - 0.004465645987498305*zp**4
    return array([[el11, el12, el13],
                  [el12, el22, el23],
                  [el13, el23, el33]], dtype=float64)/(fractional_error**2.)

def DESI_Wmat(zp):
    assert zp == 0.3

    """
    # No horizon constraint
    return array([[6991.1, 943.479, 206.799],
                  [943.479, 162.521, 22.8914],
                  [206.799, 22.8914, 6.91214]])
    """
    
    """
    # 1% constraint
    return array([[9974.89, 1944.12, 147.779],
                  [1944.12, 498.097, 3.09849],
                  [147.779, 3.09849, 8.07956]])
    """

    return array([[104227., 33552.5, -1716.54],
                  [33552.5, 11098.3, -622.12],
                  [-1716.54, -622.12, 44.9561]])

    

def get_mu(z_list, zp, deriv = None):
    r_list = n_integrate(z_list, H_inv, (0.28, -1., 0., zp), pad_list = arange(41, dtype=float64)/20.)
    if deriv == None:
        return 5.*log10((1. + z_list)*r_list)
    elif deriv == "Om":
        der_list = n_integrate(z_list, H_inv_Om, (0.28, -1., 0., zp), pad_list = arange(41, dtype=float64)/20.)
    elif deriv == "wp":
        der_list = n_integrate(z_list, H_inv_wp, (0.28, -1., 0., zp), pad_list = arange(41, dtype=float64)/20.)
    elif deriv == "wa":
        der_list = n_integrate(z_list, H_inv_wa, (0.28, -1., 0., zp), pad_list = arange(41, dtype=float64)/20.)
    return (5./log(10.))*der_list/r_list

def get_Jacobian(z_list, zp):

    jacobian = zeros([len(z_list), 4], dtype=float64)
    jacobian[:,0] = 1
    
    for i, param in enumerate(("Om", "wp", "wa")):
        jacobian[:,i+1] = get_mu(z_list, zp, deriv = param)
    return jacobian



def get_FoM(SN_Cmat, z_list, zp = 0.3, addcmb = 1, adddesi = 0, verbose = True, shift_constraint = 0.002):
    SNe_W = linalg.inv(SN_Cmat)

    jacobian = get_Jacobian(z_list, zp)

    param_W = dot(transpose(jacobian), dot(SNe_W, jacobian))

    if addcmb:
        param_W[1:,1:] += shift_parameter_Wmat(shift_constraint, zp)
    if adddesi:
        param_W[1:,1:] += DESI_Wmat(zp)


    param_C = linalg.inv(param_W)
    uncertainties = sqrt(diag(param_C))
    if verbose:
        print "Uncertainties ", uncertainties

    return 1./sqrt(linalg.det(param_C[2:4,2:4])), uncertainties

#################################################### End of FoM ###############################################








def CCM(wave, R_V): # A(wave)/A_V
    x = 1./(wave/10000.) # um^-1
    y = (x - 1.82)



    # 0.3 to 1.1
    a = (0.3 <= x)*(x <= 1.1)*0.574*(x**1.61)
    b = (0.3 <= x)*(x <= 1.1)*(-0.527*x**1.61)


    a += (1.1 < x)*(x <= 3.3)*(1. + 0.17699*y - 0.50447*y**2. - 0.02427*y**3. + 0.72085*y**4. + 0.01979*y**5. - 0.77530*y**6. + 0.32999*y**7.)
    b += (1.1 < x)*(x <= 3.3)*(1.41338*y + 2.28305*y**2. + 1.07233*y**3. - 5.38434*y**4. - 0.62251*y**5. + 5.30260*y**6. - 2.09002*y**7.)


    
    Fa = (x >= 5.9)*(-0.04473*(x - 5.9)**2. - 0.009779*(x - 5.9)**3.)
    Fb = (x >= 5.9)*(0.2130*(x - 5.9)**2. + 0.1207*(x - 5.9)**3.)

    
    a += (3.3 < x)*(x <= 8.)*(1.752 - 0.316*x - 0.104/((x - 4.67)**2. + 0.341) + Fa)
    b += (3.3 < x)*(x <= 8.)*(-3.090 + 1.825*x + 1.206/((x - 4.62)**2. + 0.263) + Fb)


    a += (8. < x)*(x <= 10.)*(-1.073 - 0.628*(x - 8.) + 0.137*(x - 8.)**2. - 0.070*(x - 8.)**3.)
    b += (8. < x)*(x <= 10.)*(13.670 + 4.257*(x - 8.) - 0.420*(x - 8.)**2. + 0.374*(x - 8.)**3.)
    

    return a + b/R_V



def write_Cmat(lines_list, sn_Cmat, FoM, flname):
    f = open(flname, 'w')
    f.write("# David Rubin's Estimated FoM: " + str(FoM) + '\n')
    f.write("#" + "\n#".join(lines_list) + '\n')
    
    for i in range(len(sn_Cmat)):
        for j in range(len(sn_Cmat)):
            f.write(str(i+1) + "  " + str(j+1) + "  " + str(sn_Cmat[i,j]) + '\n')
    f.close()

