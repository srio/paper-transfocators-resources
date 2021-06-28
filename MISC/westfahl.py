import numpy
from scipy.special import erf
import scipy.constants as codata

def tanaka_energy_spread(nu):
    num = 8 * (numpy.pi * nu)**2
    den = (2 * numpy.pi)**(3/2) * nu * erf(nu * numpy.sqrt(8 * numpy.pi**2)) + \
        numpy.exp(-8 * numpy.pi**2 * nu**2) - 1
    return num/den

def tanaka_photon_source(photon_energy, epsilon, beta, energy_spread=1e-4, Lu=2.0, period=20.0e-3, harmonic=1):
    lambda_photon = codata.h * codata.c / codata.e / photon_energy
    delta_n = 1 / (harmonic * (Lu / period))
    nu = energy_spread / delta_n
    nu_factor = (tanaka_energy_spread(energy_spread / (4 * delta_n))) ** (2/3)
    print(">>lambda_photon: ", lambda_photon)
    print(">>nu_factor: ", nu_factor)
    return numpy.sqrt(epsilon * beta + lambda_photon * Lu / (2 * numpy.pi**2) * nu_factor)

if __name__ == "__main__":

    print(tanaka_energy_spread(1e-4))

    # sirius
    beta_x = 1.5
    beta_y = 1.5
    epsilon_x = 245.0e-12
    epsilon_y = 2.45e-12
    Lu = 2.0
    period = 19e-3
    energy_spread = 1e-4


    tx = tanaka_photon_source(1000.0, epsilon_x, beta_x, energy_spread=1e-4, Lu=2.0, period=19.0e-3, harmonic=1 )
    ty = tanaka_photon_source(1000.0, epsilon_y, beta_y, energy_spread=1e-4, Lu=2.0, period=19.0e-3, harmonic=1 )

    print(1e6*tx, 1e6*ty)

