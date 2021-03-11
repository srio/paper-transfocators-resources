import numpy
import scipy.constants as codata
import xraylib


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

    element = "Be"
    density = xraylib.ElementDensity(4)

    refraction_index = xraylib.Refractive_Index(element, photon_energy/1000, density)
    refraction_index_delta = 1 - refraction_index.real
    # att_coefficient = 4*numpy.pi * (xraylib.Refractive_Index(element, photon_energy/1000, density)).imag / wave_length

    f2 = radius / 2 / refraction_index_delta
    print("which corresponds for Be to a focal length of %g m " % (f2))
    return f2, radius

def e07keV_f2_at_170m_h_calculate(fileout = "e07keV_f2_at_170m_h.dat", npoints=400,
                                  precheck=False, f1f2_grid_file=""):
    from e07keV_f2_at_170m_h import main

    if f1f2_grid_file == "":
        F1_IN = numpy.linspace(2, 100, npoints)
        F2_IN = [None] * npoints
    else:
        a = numpy.loadtxt(f1f2_grid_file)
        F1_IN = a[:,0]
        F2_IN = a[:,1]
        npoints = F1_IN.size

    F1 =   numpy.zeros(npoints)
    F2 =   numpy.zeros(npoints)
    FWHM = numpy.zeros(npoints)
    CF =   numpy.zeros(npoints)
    I0 =   numpy.zeros(npoints)

    if precheck:
        for f1 in [F1[0],F1[-1]]:
            tally = main(slit=3.62724e-05, f1=f1, f2=None)
            abscissas = tally.get_abscissas()
            spectral_density = tally.get_spectral_density_from_intensities()

            fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
            plot(abscissas, spectral_density, title="F1: %g m, FWHM: %g um" % (f1, fwhm))


    for i in range(F1.size):
        tally = main(slit=3.70441e-05, f1=F1_IN[i], f2=F2_IN[i])
        # _, occ = tally.get_occupation()
        abscissas = tally.get_abscissas()
        spectral_density = tally.get_spectral_density_from_intensities()
        fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
        # CF[i] = occ[0]
        FWHM[i] = fwhm
        I0[i] = spectral_density[spectral_density.size//2]
        F1[i] = F1_IN[i]
        if f1f2_grid_file == "":
            F2[i], _ = fit_profile(7000)
        else:
            F2[i] = F2_IN[i]
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", i, F1[i], FWHM[i], I0[i] )

    # plot(F1,FWHM)


    f = open(fileout, "w")
    f.write("#L  f1  f2  FWHM  CF  I0\n")
    for i in range(F1.size):
        f.write("%g %g %g %g %g\n" % (F1[i], F2[i], FWHM[i], CF[i], I0[i]))
    f.close()
    print("File written to disk: ", fileout)

def e07keV_f2_at_170m_v_calculate(fileout = "e07keV_f2_at_170m_h.dat", npoints=400,
                                  precheck=False, f1f2_grid_file=""):
    from e07keV_f2_at_170m_v import main

    if f1f2_grid_file == "":
        F1_IN = numpy.linspace(2, 100, npoints)
        F2_IN = [None] * npoints
    else:
        a = numpy.loadtxt(f1f2_grid_file)
        F1_IN = a[:,0]
        F2_IN = a[:,1]
        npoints = F1_IN.size

    F1 =   numpy.zeros(npoints)
    F2 =   numpy.zeros(npoints)
    FWHM = numpy.zeros(npoints)
    CF =   numpy.zeros(npoints)
    I0 =   numpy.zeros(npoints)

    if precheck:
        for i in [0,npoints-1]:
            tally = main(slit=0.000205364, f1=F1_IN[i], f2=F2_IN[i])
            abscissas = tally.get_abscissas()
            spectral_density = tally.get_spectral_density_from_intensities()

            fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
            plot(abscissas, spectral_density, title="F1: %g m, FWHM: %g um" % (F1_IN[i], fwhm))


    for i in range(F1.size):
        tally = main(slit=0.000205364, f1=F1_IN[i], f2=F2_IN[i])

        abscissas = tally.get_abscissas()
        spectral_density = tally.get_spectral_density_from_intensities()
        fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)
        print(">>>>>>>>>>>>>>>>>>> FWHM ", fwhm, FWHM.shape, F1.shape)
        # CF[i] = occ[0]
        FWHM[i] = fwhm
        I0[i] = spectral_density[spectral_density.size//2]
        F1[i] = F1_IN[i]
        if f1f2_grid_file == "":
            F2[i], _ = fit_profile(7000)
        else:
            F2[i] = F2_IN[i]
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", i, F1[i], F2[i], FWHM[i], I0[i] )



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
    plot(a[:, 0], a[:, 1], title="TRAJECTORIES", xtitle="f1 [m]", ytitle="f2 [m]")
    plot(a[:, 0], a[:, 2], title="SIZE", xtitle="f1 [m]", ytitle="size [um]")
    # plot(a[:, 0], a[:, 3], title="CF", xtitle="f1 [m]", ytitle="CF")
    # plot(a[:, 0], a[:, 4], title="I0", xtitle="f1 [m]", ytitle="I0")

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_show, set_qt
    from oasys.util.oasys_util import get_fwhm
    set_qt()


    #
    # e07keV_f2_at_170m_h_calculate(fileout="e07keV_f2_at_170m_h.dat")
    # e07keV_f2_at_170m_v_calculate(fileout="e07keV_f2_at_170m_v.dat")

    # plot_file("e07keV_f2_at_170m_h.dat")
    # plot_file("e07keV_f2_at_170m_v.dat", direction='v')

    e07keV_f2_at_170m_v_calculate(fileout="e07keV_f2_at_170m_v_f1f2_from_marco.dat",
                                  f1f2_grid_file="../GSM_MC/f1_f2_scan_7keV_f2_at_170_v.dat")

    e07keV_f2_at_170m_h_calculate(fileout="e07keV_f2_at_170m_h_f1f2_from_marco.dat",
                                  f1f2_grid_file="../GSM_MC/f1_f2_scan_7keV_f2_at_170_h.dat")








