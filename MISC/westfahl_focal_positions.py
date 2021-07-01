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

    aperture = numpy.linspace(10, 900, 100) # 200e-6
    sigma_aperture = aperture / 4.55

    p_h = 45.0
    p_v = 27.0
    f_h = 19.7
    f_v = 17.9
    q_geom_h = 1.0 / (1 / f_h - 1 / p_h)
    q_geom_v = 1.0 / (1 / f_v - 1 / p_v)

    print("position_h, position_v: ", p_h + q_geom_h, p_v + q_geom_v)
    # # esrf
    # beta_x = 6.8
    # beta_y = 2.8
    # epsilon_x = 130e-12
    # epsilon_y = 10e-12
    # Lu = 2.5
    # period = 18e-3
    # energy_spread = 1e-3
    #
    #
    #
    Photon_energy = numpy.array([3000.,9000]) #linspace(1000, 20000, 100)
    Harmonic = numpy.ones_like(Photon_energy)
    udict = {'name': "U", 'period': period*1e3, 'length': Lu, 'number': 1}


    k_array = numpy.zeros_like(Photon_energy)
    for i, ene in enumerate(Photon_energy):
        k, hh = k_value(udict, ene, electron_energy=3.0)
        Harmonic[i] = hh
        k_array[i] = k
    #
    tx =      tanaka_photon_source(Photon_energy, epsilon_x, beta_x, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )
    ty =      tanaka_photon_source(Photon_energy, epsilon_y, beta_y, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )
    tpx = tanaka_photon_divergence(Photon_energy, epsilon_x, beta_x, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )
    tpy = tanaka_photon_divergence(Photon_energy, epsilon_y, beta_y, energy_spread=energy_spread, Lu=Lu, period=period, harmonic=Harmonic )


    print("Sizes [um] H at 3, 8 keV: ", 1e6 * tx)
    print("Sizes [um] V at 3, 8 keV: ", 1e6 * ty)
    #
    #
    # #
    # # elleaume
    # #
    #
    # sx = numpy.sqrt(epsilon_x * beta_x)
    # sy = numpy.sqrt(epsilon_y * beta_y)
    # spx = numpy.sqrt(epsilon_x / beta_x)
    # spy = numpy.sqrt(epsilon_y / beta_y)
    #
    # # for i, ene in enumerate(Photon_energy):
    # #     # k, hh = k_value(u18, ene, electron_energy=3.0)
    # #     # Harmonic[i] = hh
    sr = tanaka_photon_source(Photon_energy, 0, 1, energy_spread=energy_spread, Lu=2.0, period=19.0e-3,
                              harmonic=Harmonic)
    srp = tanaka_photon_divergence(Photon_energy, 0, 1, energy_spread=energy_spread, Lu=2.0, period=19.0e-3,
                              harmonic=Harmonic)


    CF_h = sr * srp / (tx * tpx)
    CF_v = sr * srp / (ty * tpy)
    beta_GSM_h = CF_h / numpy.sqrt(1 - CF_h)
    beta_GSM_v = CF_v / numpy.sqrt(1 - CF_v)

    print("CF H,V: ", CF_h, CF_v)
    print("beta GSM H,V: ", beta_GSM_h, beta_GSM_v)

    sigma_GSM = tx
    xi_GSM = beta_GSM_h * sigma_GSM
    print("xi GSM 3, 9 keV: ", xi_GSM)

    wavelength = codata.h * codata.c / codata.e / Photon_energy
    Z_R = numpy.pi * sigma_GSM**2 / wavelength * (1 + (sigma_GSM/xi_GSM)**2 )**(-1/2)

    print("Z_R 3, 9 keV: ", Z_R)

    first_term_num = p_h / f_h - 1
    first_term_den = (p_h / f_h - 1)**2 + (Z_R / f_h)**2
    q = f_h * (1 + first_term_num / first_term_den)

    print("p_h + q: ", p_h + q)

    #
    #
    #
    Aperture = numpy.linspace(50e-6, 900e-6, 100)
    Sigma_aperture = Aperture / 4.55 * numpy.ones_like(Aperture)

    sigmap_GSM = wavelength / 4 / numpy.pi / sigma_GSM
    print(">>> sigma_GSM, sigmap_GSM at 3, 9 keV", sigma_GSM, sigmap_GSM)
    Sigma_beam_before_aperture_3 = numpy.sqrt( (sigmap_GSM[0] * p_h) ** 2 + sigmap_GSM[0] ** 2) * numpy.ones_like(Aperture)
    Sigma_beam_before_aperture_9 = numpy.sqrt( (sigmap_GSM[1] * p_h) ** 2 + sigmap_GSM[1] ** 2) * numpy.ones_like(Aperture)

    Sigma_beam_after_aperture_3 = numpy.sqrt(1.0 / (1 / Sigma_beam_before_aperture_3**2 + 1 / Sigma_aperture**2))
    Sigma_beam_after_aperture_9 = numpy.sqrt(1.0 / (1 / Sigma_beam_before_aperture_9**2 + 1 / Sigma_aperture**2))


    print("shapes: ", Sigma_beam_after_aperture_3.shape, Sigma_beam_after_aperture_9.shape)

    beta_GSM_h_propagated_3 = beta_GSM_h[0] * Sigma_beam_before_aperture_3 / Sigma_beam_after_aperture_3
    beta_GSM_h_propagated_9 = beta_GSM_h[1] * Sigma_beam_before_aperture_9 / Sigma_beam_after_aperture_9

    Z_R_aperture_3 = numpy.pi * Sigma_beam_after_aperture_3 ** 2 / wavelength[0] # * (1 + (1.0 / beta_GSM_h_propagated_3) ** 2 )**(-1/2)
    Z_R_aperture_9 = numpy.pi * Sigma_beam_after_aperture_9 ** 2 / wavelength[1] # * (1 + (1.0 / beta_GSM_h_propagated_9) ** 2) **(-1/2)

    # Z_R_aperture_3 = numpy.zeros_like(Aperture)
    # Z_R_aperture_9 = numpy.zeros_like(Aperture)

    print("Wavelength: ", wavelength)
    print("Z_R_aperture: ", Z_R_aperture_3, Z_R_aperture_9)


    ratio_Z_R_over_f_h_3 = Z_R_aperture_3 / f_h
    ratio_Z_R_over_f_h_9 = Z_R_aperture_9 / f_h
    ratio_Z_R_over_f_h_G = numpy.zeros_like(Aperture)


    # ratio_Z_R_over_f_h_3 = 0.1 * numpy.ones_like(Aperture)
    # ratio_Z_R_over_f_h_9 = 1.0 * numpy.ones_like(Aperture)


    first_term_den_3 = (p_h / f_h - 1)**2 + (ratio_Z_R_over_f_h_3)**2
    Q_3 = f_h * (1 + first_term_num / first_term_den_3)
    first_term_den_9 = (p_h / f_h - 1)**2 + (ratio_Z_R_over_f_h_9)**2
    Q_9 = f_h * (1 + first_term_num / first_term_den_9)
    first_term_den_G = (p_h / f_h - 1)**2 + (ratio_Z_R_over_f_h_G)**2
    Q_G = f_h * (1 + first_term_num / first_term_den_G)

    from srxraylib.plot.gol import plot
    print(Q_3)
    plot(1e6 * Aperture, p_h + Q_3,
         1e6 * Aperture, p_h + Q_9,
         1e6 * Aperture, p_h + Q_G,
         1e6 * Aperture, p_h + f_h * numpy.zeros_like(Aperture),
         legend=["3 keV", "9 keV", "geometric", "at slit"])

    plot(1e6 * Aperture, 1e6 * Sigma_beam_before_aperture_3,
         1e6 * Aperture, 1e6 * Sigma_beam_after_aperture_3,
         1e6 * Aperture, 1e6 * Sigma_aperture,
        legend=["before", "after", "aperture"])

    plot(
         # 1e6 * Aperture, Z_R_aperture_3,
         # 1e6 * Aperture, Z_R_aperture_9,
         1e6 * Aperture, ratio_Z_R_over_f_h_3,
         1e6 * Aperture, ratio_Z_R_over_f_h_9,
        legend=["Z_R/f 3", "Z_R/f 9"])

    plot(
         # 1e6 * Aperture, Z_R_aperture_3,
         # 1e6 * Aperture, Z_R_aperture_9,
         1e6 * Aperture, 1 / first_term_den_3,
         1e6 * Aperture, 1 / first_term_den_9,
        legend=["1/den 3", "1/den 9"])

    # # srp /= 2.0 #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #
    # Sx = numpy.sqrt( sr ** 2  + sx ** 2 )
    # Sy = numpy.sqrt( sr ** 2  + sy ** 2 )
    # Spx = numpy.sqrt(srp ** 2 + spx ** 2)
    # Spy = numpy.sqrt(srp ** 2 + spy ** 2)
    # print(sx, sy)
    #
    # from srxraylib.plot.gol import plot
    # plot(Photon_energy, 1e6 * tx,
    #      Photon_energy, 1e6 * ty,
    #      Photon_energy, 1e6 * Sx,
    #      Photon_energy, 1e6 * Sy,
    #      # Photon_energy, 1e6 * Sx,
    #      # Photon_energy, 1e6 * Sy,
    #      yrange=[0,35],
    #      legend=["H Tanaka", "V Tanaka", "H Elleaume", "V Elleaume"],
    #      title="Sizes", show=0)
    #
    #
    # plot(Photon_energy, 1e6 * tpx,
    #      Photon_energy, 1e6 * tpy,
    #      Photon_energy, 1e6 * Spx,
    #      Photon_energy, 1e6 * Spy,
    #      yrange=[0, 25],
    #      legend=["H Tanaka","V Tanaka","H Elleaume","V Elleaume"],
    #      title="Divergences", show=1)



