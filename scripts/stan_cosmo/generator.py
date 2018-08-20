import numpy as np


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

