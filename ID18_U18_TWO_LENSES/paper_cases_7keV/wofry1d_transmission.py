import numpy
import h5py
from srxraylib.plot.gol import plot_image, plot, set_qt
set_qt()

do_plot=0

def get_integral(file, do_plot=do_plot):

    f = h5py.File(file, 'r')

    image = f["/entry/data0/image"][()]
    x = f["/entry/data0/x"][()]
    y = f["/entry/data0/y"][()]

    f.close()


    print(image.shape, x.shape, y.shape)


    if do_plot:
        plot_image(image.T, x, y)


    integral = image.sum() * (x[1] - x[0]) * (y[1] - y[0])



    return integral

def get_integral_from_1d_profiles(fileh, filev, do_plot=do_plot):

    ah = numpy.loadtxt(fileh)
    av = numpy.loadtxt(filev)

    x = ah[:,0]
    y = av[:,0]

    if do_plot:
        plot(ah[:,0],ah[:,1],
             av[:,0],av[:,1],
             legend=['h','v'])

    image = numpy.outer(ah[:,1], av[:,1])

    integral = image.sum() * (x[1] - x[0]) * (y[1] - y[0])

    return integral

def get_integral_from_1d_wavefronts(wfr1, wfr2, do_plot=do_plot):

    f = h5py.File(wfr1, 'r')
    intensity1 = f["/wfr/intensity/wfr_intensity"][()]
    x          = f["/wfr/intensity/axis_x"][()]
    f.close()


    f = h5py.File(wfr2, 'r')
    intensity2 = f["/wfr/intensity/wfr_intensity"][()]
    y          = f["/wfr/intensity/axis_x"][()]
    f.close()

    if do_plot:
        plot(x, intensity1,
             y, intensity2,
             legend=['h','v'])

    image = numpy.outer(intensity1, intensity2)

    integral = image.sum() * (x[1] - x[0]) * (y[1] - y[0])

    return integral

if __name__ == "__main__":

    # #
    # # OLD (writing by hand the 2d profile from 1d->2d widget)
    # #
    #
    # dir = "/users/srio/Oasys/"
    # integral0 = get_integral(dir+"case1_source.h5")
    # print(integral0)
    # integral1 = get_integral(dir+"case1_slit.h5")
    # print(integral1)
    # integral2 = get_integral(dir+"case1_lens1.h5")
    # print(integral2)
    # integral3 = get_integral(dir+"case1_lens2.h5")
    # print(integral3)
    #
    # print("Absorption slit, lens1, lens2 : %3.1f, %3.1f, %3.1f " % (100 * (integral0 - integral1) / integral0,
    #                                                   100 * (integral1 - integral2) / integral1,
    #                                                   100 * (integral2 - integral3) / integral2,) )

    #
    # running CMD script with 50 modes
    #
    dir = "/users/srio/Oasys/"
    case = 4

    integral0 = get_integral_from_1d_profiles(dir + "case%dh_wofry_source_spectral_density.dat" % case,
                                              dir + "case%dv_wofry_source_spectral_density.dat" % case)
    # print(integral0)

    integral1 = get_integral_from_1d_profiles(dir + "case%dh_wofry_slit_spectral_density.dat" % case,
                                              dir + "case%dv_wofry_slit_spectral_density.dat" % case)
    # print(integral1)

    integral2 = get_integral_from_1d_profiles(dir + "case%dh_wofry_lens1_spectral_density.dat" % case,
                                              dir + "case%dv_wofry_lens1_spectral_density.dat" % case)
    # print(integral2)

    integral3 = get_integral_from_1d_profiles(dir + "case%dh_wofry_lens2_spectral_density.dat" % case,
                                              dir + "case%dv_wofry_lens2_spectral_density.dat" % case)
    # print(integral3)


    print("CASE %d PARTIAL COH Absorption slit, lens1, lens2 : %3.1f, %3.1f, %3.1f " % (case, 100 * (integral0 - integral1) / integral0,
                                                      100 * (integral1 - integral2) / integral1,
                                                      100 * (integral2 - integral3) / integral2,) )

    # print("Absorption slit %3.1f " % (100 * (integral0 - integral1) / integral0 ))


    #
    # from wavefronts
    #

    dir = "/users/srio/Oasys/"

    integral0 = get_integral_from_1d_wavefronts(dir + "case%dh_wofry_source.h5" % case,
                                                dir + "case%dv_wofry_source.h5" % case)
    # print(integral0)

    integral1 = get_integral_from_1d_wavefronts(dir + "case%dh_wofry_slit.h5" % case,
                                                dir + "case%dv_wofry_slit.h5" % case)
    # print(integral1)

    integral2 = get_integral_from_1d_wavefronts(dir + "case%dh_wofry_lens1.h5" % case,
                                                dir + "case%dv_wofry_lens1.h5" % case)
    # print(integral2)

    integral3 = get_integral_from_1d_wavefronts(dir + "case%dh_wofry_lens2.h5" % case,
                                                dir + "case%dv_wofry_lens2.h5" % case)
    # print(integral3)


    print("CASE %d FULL COH Absorption slit, lens1, lens2 : %3.1f, %3.1f, %3.1f " % (case,100 * (integral0 - integral1) / integral0,
                                                      100 * (integral1 - integral2) / integral1,
                                                      100 * (integral2 - integral3) / integral2,) )

    # print("Absorption slit %3.1f " % (100 * (integral0 - integral1) / integral0 ))







