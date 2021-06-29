# source sizes as in https://doi.org/10.1107/S1600577517003058

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
    delta_n = 1.0 / (harmonic * (Lu // period))
    nu_factor = (tanaka_energy_spread(energy_spread / (4 * delta_n))) ** (2/3)
    # print(">>lambda_photon: ", lambda_photon)
    # print(">>nu_factor**(3/2): ", nu_factor)
    return numpy.sqrt(epsilon * beta + lambda_photon * Lu / (2 * numpy.pi**2) * nu_factor)

def tanaka_photon_divergence(photon_energy, epsilon, beta, energy_spread=1e-4, Lu=2.0, period=20.0e-3, harmonic=1):
    print(Harmonic)
    lambda_photon = codata.h * codata.c / codata.e / photon_energy
    delta_n = 1.0 / (harmonic * (Lu // period))
    nu = energy_spread / delta_n
    nu_factor = (tanaka_energy_spread(energy_spread / delta_n))
    # print(">>lambda_photon: ", lambda_photon)
    # print(">>nu_factor: ", nu_factor)
    return numpy.sqrt(epsilon / beta + lambda_photon / (2 * Lu ) * nu_factor)

# def k_value(und_dict, photon_energy, electron_energy=6.0):
#
#     gamma = electron_energy* 1e9 / (codata.m_e *  codata.c**2 / codata.e)
#     lambdan = codata.h*codata.c/codata.e*1e10 / photon_energy # in A
#
#     harm_number = -1
#     KK = -1
#     while KK < 0:
#         harm_number += 2
#         lambda1 = lambdan * harm_number
#         KK = ( 2*( (1e-10*lambda1)/(1e-3*und_dict['period']) *2*gamma*gamma  - 1))
#
#     return numpy.sqrt(KK),harm_number

def k_value(und_dict, photon_energy, electron_energy=6.0, kmin=0.0, kmax=2):

        gamma = electron_energy * 1e9 / (codata.m_e * codata.c ** 2 / codata.e)
        lambdan = codata.h * codata.c / codata.e * 1e10 / photon_energy  # in A

        harm_number = -1
        KK = -1
        igood = 0

        while igood == 0:
            harm_number += 2
            lambda1 = lambdan * harm_number
            KK = (2 * ((1e-10 * lambda1) / (1e-3 * und_dict['period']) * 2 * gamma * gamma - 1))
            if KK < 0:
                igood = 0
            else:
                print(">>>>>>>>>>>>>>>>>>>>>>>>", numpy.sqrt(KK))
                kv = numpy.sqrt(KK)
                if  kv <= kmax and kv >= kmin:
                    igood = 1
                else:
                    if kv > kmax:
                        igood = 1
                        kv = 0
                        harm_number = 0
                    else:
                        igood = 0

        return kv, harm_number

def get_sigmas_radiation(photon_energy,undulator_length):
    lambdan = 1e-10 * codata.h*codata.c/codata.e*1e10 / photon_energy # in m
    return 2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),0.69*numpy.sqrt(lambdan/undulator_length),

if __name__ == "__main__":

    # sirius
    beta_x = 1.5
    beta_y = 1.5
    epsilon_x = 245.0e-12
    epsilon_y = 2.45e-12
    Lu = 2.0
    period = 19e-3
    energy_spread = 1e-3

    # # esrf
    # beta_x = 6.8
    # beta_y = 2.8
    # epsilon_x = 130e-12
    # epsilon_y = 10e-12
    # Lu = 2.5
    # period = 18e-3
    # energy_spread = 1e-3



    Photon_energy = numpy.linspace(1000, 20000, 100)
    Harmonic = numpy.ones_like(Photon_energy)
    udict = {'name': "U", 'period': period*1e3, 'length': Lu, 'number': 1}


    k_array = numpy.zeros_like(Photon_energy)
    for i, ene in enumerate(Photon_energy):
        k, hh = k_value(udict, ene, electron_energy=3.0)
        Harmonic[i] = hh
        k_array[i] = k

    tx =      tanaka_photon_source(Photon_energy, epsilon_x, beta_x, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )
    ty =      tanaka_photon_source(Photon_energy, epsilon_y, beta_y, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )
    tpx = tanaka_photon_divergence(Photon_energy, epsilon_x, beta_x, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )
    tpy = tanaka_photon_divergence(Photon_energy, epsilon_y, beta_y, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )


    #
    # elleaume
    #

    sx = numpy.sqrt(epsilon_x * beta_x)
    sy = numpy.sqrt(epsilon_y * beta_y)
    spx = numpy.sqrt(epsilon_x / beta_x)
    spy = numpy.sqrt(epsilon_y / beta_y)

    # for i, ene in enumerate(Photon_energy):
    #     # k, hh = k_value(u18, ene, electron_energy=3.0)
    #     # Harmonic[i] = hh
    sr, srp = get_sigmas_radiation(Photon_energy, Lu)
    # srp /= 2.0 #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    Sx = numpy.sqrt( sr ** 2  + sx ** 2 )
    Sy = numpy.sqrt( sr ** 2  + sy ** 2 )
    Spx = numpy.sqrt(srp ** 2 + spx ** 2)
    Spy = numpy.sqrt(srp ** 2 + spy ** 2)
    print(sx, sy)

    from srxraylib.plot.gol import plot
    plot(Photon_energy, 1e6 * tx,
         Photon_energy, 1e6 * ty,
         Photon_energy, 1e6 * Sx,
         Photon_energy, 1e6 * Sy,
         # Photon_energy, 1e6 * Sx,
         # Photon_energy, 1e6 * Sy,
         yrange=[0,35],
         legend=["H Tanaka", "V Tanaka", "H Elleaume", "V Elleaume"],
         title="Sizes", show=0)


    plot(Photon_energy, 1e6 * tpx,
         Photon_energy, 1e6 * tpy,
         Photon_energy, 1e6 * Spx,
         Photon_energy, 1e6 * Spy,
         yrange=[0, 25],
         legend=["H Tanaka","V Tanaka","H Elleaume","V Elleaume"],
         title="Divergences", show=1)



