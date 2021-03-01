import numpy

def parse_data_marco(display_original_fh=False):
    directions = ['v','h']  # do not change order!
    focii = ['large', 'small', ]
    energies = [7, 10, 15, 20, 35]
    distances = [170, 192]

    import datastorage
    data = datastorage.read("summary_of_GSM_results_for_Manuel.npz")


    for distance in distances:
        for energy in energies:
            for focus in focii:
                label = "e%02dkeV_f2_at_%dm_%s" % (energy, distance, focus)

                tmp = data[label]['v']
                v_slit      = tmp["p035m_hard_aperture"]
                v_f1        = tmp["p066m_f1_used"]
                v_f2        = tmp["pf2pos_f2_used"]
                v_f1d       = tmp["p066m_f1_desired"]
                v_f2d       = tmp["pf2pos_f2_desired"]
                v_size200   = tmp["p200m_fwhm_beam"]
                v_sizewaist = tmp["waist_fwhm_size"]
                v_poswaist  = tmp["waist_pos"]


                tmp = data[label]['h']
                h_slit = tmp["p035m_hard_aperture"]
                h_f1 = tmp["p066m_f1_used"]
                h_f2 = tmp["pf2pos_f2_used"]
                h_f1d = tmp["p066m_f1_desired"]
                h_f2d = tmp["pf2pos_f2_desired"]
                h_size200 = tmp["p200m_fwhm_beam"]
                h_sizewaist = tmp["waist_fwhm_size"]
                h_poswaist = tmp["waist_pos"]

                if display_original_fh:
                    print(">>>>%s %s :  slit: %g m, f1: %g (%g), f2: %g (%g), size: %g um (%g at %g m)" %
                          (label, 'h', h_slit, h_f1, h_f1d, h_f2, h_f2d, 1e6*h_size200, 1e6*h_sizewaist, h_poswaist))
                else:
                    print(">>>>%s %s :  slit: %g m, f1: %g (%g), f2: %g (%g), size: %g um (%g at %g m)" %
                          (label, 'h', h_slit,
                           1/(1/h_f1+1/v_f1), 1/(1/h_f1d+1/v_f1d), 1/(1/h_f2+1/v_f2), 1/(1/h_f2d+1/v_f2d),
                           1e6*h_size200, 1e6*h_sizewaist, h_poswaist))


                print(">>>>%s %s :  slit: %g m, f1: %g (%g), f2: %g (%g), size: %g um (%g at %g m)" %
                      (label, 'v', v_slit, v_f1, v_f1d, v_f2, v_f2d, 1e6 * v_size200, 1e6 * v_sizewaist, v_poswaist))




def get_data_marco(label="e07keV_f2_at_170m",direction='h'):
    # directions = ['h','v']
    # focii = ['large', 'small', ]
    # energies = [7, 10, 15, 20, 35]
    # distances = [170, 192]

    import datastorage
    data = datastorage.read("summary_of_GSM_results_for_Manuel.npz")

    # # print("%s  %s" % (label, direction))
    # tmp = data[label][direction]
    # for key in tmp.keys():
    #     print(key, tmp[key])
    #     slit      = tmp["p035m_hard_aperture"]
    #     f1        = tmp["p066m_f1_used"]
    #     f2        = tmp["pf2pos_f2_used"]
    #     f1w         = tmp["p066m_f1_desired"]
    #     f2w        = tmp["pf2pos_f2_desired"]
    #     size200   = tmp["p200m_fwhm_beam"]
    #     sizewaist = tmp["waist_fwhm_size"]
    #     poswaist  = tmp["waist_pos"]
    # # print(">>>>%s %s :  slit: %g m, f1: %g, f2: %g, size: %g um (%g at %g m)" %
    # #       (label, direction, slit, f1, f2, 1e6*size200, 1e6*sizewaist, poswaist))

    tmp = data[label]['v']
    v_slit = tmp["p035m_hard_aperture"]
    v_f1 = tmp["p066m_f1_used"]
    v_f2 = tmp["pf2pos_f2_used"]
    v_f1d = tmp["p066m_f1_desired"]
    v_f2d = tmp["pf2pos_f2_desired"]
    v_size200 = tmp["p200m_fwhm_beam"]
    v_sizewaist = tmp["waist_fwhm_size"]
    v_poswaist = tmp["waist_pos"]

    tmp = data[label]['h']
    h_slit = tmp["p035m_hard_aperture"]
    h_f1 = tmp["p066m_f1_used"]
    h_f2 = tmp["pf2pos_f2_used"]
    h_f1d = tmp["p066m_f1_desired"]
    h_f2d = tmp["pf2pos_f2_desired"]
    h_size200 = tmp["p200m_fwhm_beam"]
    h_sizewaist = tmp["waist_fwhm_size"]
    h_poswaist = tmp["waist_pos"]

    if direction == 'h':
        return h_slit, \
               1 / (1 / h_f1 + 1 / v_f1), 1 / (1 / h_f2 + 1 / v_f2), \
               1 / (1 / h_f1d + 1 / v_f1d), 1 / (1 / h_f2d + 1 / v_f2d), \
               h_size200, h_sizewaist # , h_poswaist
    else:
        return v_slit, v_f1, v_f2, v_f1d, v_f2d, v_size200, v_sizewaist #, v_poswaist
    #
    #
    # return slit, f1, f2, f1w, f2w, size200, sizewaist

def fit_profile(photon_energy):
    a = numpy.loadtxt('profile1D.dat')
    n = a.shape[0]
    w = n // 20

    x = a[(n // 2 - w):(n // 2 + w), 0] * 1e-6
    y = a[(n // 2 - w):(n // 2 + w), 1] * 1e-6

    yder = numpy.gradient(y, x)
    coeff = numpy.polyfit(x, yder, 1)

    print("\n\n\n ==========  fitted radius in the profile center : ")
    radius = 2e6 / coeff[0]
    print("fitted lens (with two curved sides) of radius = %g m " % (radius))


    import scipy.constants as codata
    import xraylib


    element = "Be"
    density = xraylib.ElementDensity(4)

    refraction_index = xraylib.Refractive_Index(element, photon_energy/1000, density)
    refraction_index_delta = 1 - refraction_index.real
    # att_coefficient = 4*numpy.pi * (xraylib.Refractive_Index(element, photon_energy/1000, density)).imag / wave_length

    f2 = radius / 2 / refraction_index_delta
    print("which corresponds for Be to a focal length of %g m " % (f2))
    return f2

def e07keV_f2_at_170m_h_calculate(fileout = "e07keV_f2_at_170m_h.dat"):
    from e07keV_f2_at_170m_h import main

    # F1 = numpy.array([30.53673,102.6058])
    F1 = numpy.linspace(7.0,500.0,500)
    F2 = numpy.zeros_like(F1)
    FWHM = numpy.zeros_like(F1)
    CF = numpy.zeros_like(F1)
    I0 = numpy.zeros_like(F1)

    for f1 in [F1[0],F1[-1]]:
        tally = main(slit=3.62724e-05, f1=f1, f2=None)
        # tally.plot_spectral_density()
        # _, cf = tally.get_occupation()
        abscissas = tally.get_abscissas()
        spectral_density = tally.get_spectral_density_from_intensities()

        fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
        plot(abscissas, spectral_density, title="F1: %g m, FWHM: %g um" % (f1, fwhm))


    for i in range(F1.size):
        tally = main(slit=3.62724e-05, f1=F1[i], f2=None)
        # _, occ = tally.get_occupation()
        abscissas = tally.get_abscissas()
        spectral_density = tally.get_spectral_density_from_intensities()
        fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
        # CF[i] = occ[0]
        FWHM[i] = fwhm
        I0[i] = spectral_density[spectral_density.size//2]
        F2[i] = fit_profile(7000)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", i, F1[i], FWHM[i], I0[i] )




    plot(F1,FWHM)


    f = open(fileout, "w")
    f.write("#L  f1  f2  FWHM  CF  I0\n")
    for i in range(F1.size):
        f.write("%g %g %g %g %g\n" % (F1[i], F2[i], FWHM[i], CF[i], I0[i]))
    f.close()
    print("File written to disk: ", fileout)

def e07keV_f2_at_170m_v_calculate(fileout = "e07keV_f2_at_170m_h.dat"):
    from e07keV_f2_at_170m_v import main

    F1 = numpy.array([47.8825,1000.0])
    F1 = numpy.concatenate((numpy.linspace(15.0,500.0,500),numpy.linspace(501,1500.0,5)))
    F2 = numpy.zeros_like(F1)
    FWHM = numpy.zeros_like(F1)
    CF = numpy.zeros_like(F1)
    I0 = numpy.zeros_like(F1)

    for f1 in [F1[0],F1[-1]]:
        tally = main(slit=201e-6, f1=f1, f2=None)
        # tally.plot_spectral_density()
        # _, cf = tally.get_occupation()
        abscissas = tally.get_abscissas()
        spectral_density = tally.get_spectral_density_from_intensities()

        fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
        plot(abscissas, spectral_density, title="F1: %g m, FWHM: %g um" % (f1, fwhm))


    for i in range(F1.size):
        tally = main(slit=201e-6, f1=F1[i])
        # _, occ = tally.get_occupation()
        abscissas = tally.get_abscissas()
        spectral_density = tally.get_spectral_density_from_intensities()
        fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
        # CF[i] = occ[0]
        FWHM[i] = fwhm
        I0[i] = spectral_density[spectral_density.size//2]
        F2[i] = fit_profile(7000)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", i, F1[i], FWHM[i], I0[i] )




    plot(F1,FWHM)


    f = open(fileout, "w")
    f.write("#L  f1  f2  FWHM  CF  I0\n")
    for i in range(F1.size):
        f.write("%g %g %g %g %g\n" % (F1[i], F2[i], FWHM[i], CF[i], I0[i]))
    f.close()
    print("File written to disk: ", fileout)

def plot_file(fileout, direction='h'):
    # fileout = "e07keV_f2_at_170m_h.dat"
    a = numpy.loadtxt(fileout, skiprows=1)
    print(a.shape)


    #
    # plot(a[:, 0], a[:, 3], title="CF", xtitle="f1 [m]", ytitle="CF", xrange=[0,100])
    # plot(a[:, 0], a[:, 4], title="I0", xtitle="f1 [m]", ytitle="I0", xrange=[0,100])


    slitl, f1l, f2l, f1wl, f2wl, size200l, sizewaistl = get_data_marco("e07keV_f2_at_170m_large", direction=direction)
    slits, f1s, f2s, f1ws, f2ws, size200s, sizewaists = get_data_marco("e07keV_f2_at_170m_small", direction=direction)
    print(">>>>>>>Small", slits, f1s, f2s, f1ws, f2ws, size200s, sizewaists)
    print(">>>>>>>Large", slitl, f1l, f2l, f1wl, f2wl, size200l, sizewaistl)

    if f1s > 1000: f1s = 1000
    if f1ws > 1000: f1ws = 1000
    # if F1l > 1000: F1s = 1000

    if direction == 'h':
        xrange = [0, 115]
        yrange = [15, 35]
    else:
        xrange = [0, 1015]
        yrange = [15, 55]

    f = plot(a[:, 0], a[:, 1],
             [f1l, f1s], [f2l, f2s],
             [f1wl, f1ws], [f2wl, f2ws],
             title=fileout + " trajectories", xtitle="f1 [m]", ytitle="f2 [m]",
             xrange=xrange, yrange=yrange,
             show=1,
             marker=["None",'+','o'],
             legend=["Wofry","Marco (TF choice)","Marco (ideal choice)"])


    g = plot(a[:, 0], a[:, 2],
             [f1l, f1s], [size200l * 1e6, size200s * 1e6],
             [f1wl, f1ws], [sizewaistl * 1e6, sizewaists * 1e6],
             title=fileout + " FWHM", xtitle="f1 [m]", ytitle="FWHM [um]",
             xrange=xrange,
             # xrange=[0,1015],
             show=0,
             marker=["None", '+', 'o'],
             legend=["Wofry", "Marco (TF choice @200m)","Marco (TF choice @waist)"])
    g[1].grid()
    plot_show()


if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_show, set_qt
    from oasys.util.oasys_util import get_fwhm
    set_qt()

    parse_data_marco()

    # e07keV_f2_at_170m_h_calculate(fileout="e07keV_f2_at_170m_h.dat")
    plot_file("e07keV_f2_at_170m_h.dat")


    # e07keV_f2_at_170m_v_calculate(fileout="e07keV_f2_at_170m_v.dat")
    plot_file("e07keV_f2_at_170m_v.dat", direction='v')










