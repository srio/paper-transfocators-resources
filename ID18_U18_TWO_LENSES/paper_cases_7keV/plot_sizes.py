import numpy
from srxraylib.plot.gol import plot, set_qt
from silx.io.specfile import SpecFile
import xraylib
import matplotlib.pylab as plt
from scipy.signal import savgol_filter


def get_f2(f1=28.2,
           position_source=0.0,
           position_lens1=65.0,
           position_lens2= 170.0,
           position_sample=200.0,
           verbose=True):

    p1 = position_lens1 - position_source
    q1 = 1 / (1 / f1 - 1 / p1)


    p2 = position_lens2 - (p1 + q1)
    q2 = position_sample - position_lens2

    f2 = 1.0 / (1 / p2 + 1 / q2)

    if verbose:
        D = position_lens2 - position_lens1
        print("D: %g, q1+p2: %g" % (D, q1+p2))
        print("p1: %g" % p1)
        print("q1: %g" % q1)
        print("p2: %g" % p2)
        print("q2: %g" % q2)
        print("D: %g, Sum: %g" % (D, p1+q1+p2+q2))

    M = (q1 / q1) * (q2 / p2)
    return f2, M

def plot_one_case(direction='h',ii=0):


    if direction == 'h':
        APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]
        NFILES = [200, 200, 200, 200]
        # ii = 2
    else:
        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]
        NFILES = [200, 200, 200, 200]
        # ii = 1

    aperture = APERTURE[ii]
    nfiles = NFILES[ii]


    #
    #
    #

    # for aperture in APERTURE:
    if True:

        Index = numpy.arange(nfiles)
        if direction == 'h':
            sourcesize = 70.57e-6
        else:
            sourcesize = 15.02e-6

        a = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (aperture*1e6, direction))
        F1 = a[0:nfiles, 0].copy()
        F2 = a[0:nfiles, 1].copy()


        #
        # read files with sizes
        #
        FWHM = []
        for index in Index:
            try:
                filename = "sizes_slit%g_%s/case7keV_%s_spectral_density_%03d.dat" % (aperture*1e6, direction, direction, index)
                # filename = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/sizes_slit40.3_h/case7keV_h_spectral_density_000.dat"
                print(">>>>>>>>> opening file: ",filename)
                sf = SpecFile(filename)
                s1 = sf[0]
                fwhm = s1.scan_header_dict["UFWHM"]
                FWHM.append(float(fwhm))
            except:
                FWHM.append(0)


        F2theory1 = []
        F2theory2 = []
        Msource_at_id = []
        Msource_at_slit = []

        for index in Index:
            ff_source_at_id, mm_source_at_id = get_f2(f1=F1[index],
                                     position_source=0.0,  # source at source
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=False)

            ff_source_at_slit, mm_source_at_slit = get_f2(f1=F1[index],
                                     position_source=36.0, #0.0,  # source at slit
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=False)

            Msource_at_id.append(mm_source_at_id)
            Msource_at_slit.append(mm_source_at_slit)

            F2theory1.append(ff_source_at_id)
            F2theory2.append(ff_source_at_slit)


        print("F:", len(F1), len(F2))


        F2theory1smooth = numpy.array(F2theory1)
        F2theory2smooth = numpy.array(F2theory2)
        Msource_at_id = numpy.array(Msource_at_id)
        Msource_at_slit = numpy.array(Msource_at_slit)


        fig, ax = plot(numpy.array(F1), numpy.array(F2),
            numpy.array(F1), F2theory1smooth,
            numpy.array(F1), F2theory2smooth,
            legend=["Wofry1D","Geometrical Optics smoothed (source at ID)","Geometrical Optics (source at slit)"],
            marker=[None,None,None],
            linestyle=['-',':',':'],
            yrange=[10,60],
            xtitle="F1 [m]", ytitle="F2 [m]", title="%s trajectories (slit=%g um)" % (direction, aperture*1e6),
            show=0)

        ax.xaxis.grid()
        ax.yaxis.grid()


        fig2, ax2 = plot(numpy.array(F1), numpy.array(FWHM),
             numpy.array(F1), Msource_at_id * sourcesize * 1e6,
             numpy.array(F1), Msource_at_slit * aperture * 1e6,
             marker=[None,'.','.'],
             linestyle=['-',"None","None"],
             yrange=[0,numpy.array(FWHM).max()*1.1],
             legend=["Wofry1D","Geometrical optics (source at ID)","Geometrical optics (source at slit)"],
             xtitle="F1 [m]", ytitle="FWHM [um]", title="%s Sizes (slit=%g um)" % (direction, aperture*1e6),
             show=0)

        ax2.xaxis.grid()
        ax2.yaxis.grid()

        plt.show()

        # if filename is not None:
        #     plt.savefig(filename)
        #     print("File written to disk: %s" % filename)


def plot_all_cases(direction='h',graphic_filename=None):


    if direction == 'h':
        APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]
        NFILES = [200, 200, 200, 200]
        # ii = 2
    else:
        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]
        NFILES = [200, 200, 200, 200]
        # ii = 1

    # aperture = APERTURE[ii]
    # nfiles = NFILES[ii]




    #
    #
    #

    # for aperture in APERTURE:
    if True:

        Index = numpy.arange(200)
        if direction == 'h':
            sourcesize = 70.57e-6
            yrange=[1,500] # [0,75]
        else:
            sourcesize = 15.02e-6
            yrange = [1,500] # [0, 200]

        # a = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (aperture*1e6, direction))

        a_0 = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[0], direction))
        a_1 = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[1], direction))
        a_2 = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[2], direction))
        a_3 = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[3], direction))


        F1_0 = a_0[:, 0].copy()
        F2_0 = a_0[:, 1].copy()

        F1_1 = a_1[:, 0].copy()
        F2_1 = a_1[:, 1].copy()

        F1_2 = a_2[:, 0].copy()
        F2_2 = a_2[:, 1].copy()

        F1_3 = a_3[:, 0].copy()
        F2_3 = a_3[:, 1].copy()

        #
        # read files with sizes
        #
        FWHM_0 = []

        for index in Index:
            try:
                filename = "sizes_slit%g_%s/case7keV_%s_spectral_density_%03d.dat" % (APERTURE[0]*1e6, direction, direction, index)
                # filename = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/sizes_slit40.3_h/case7keV_h_spectral_density_000.dat"
                print(">>>>>>>>> opening file: ",filename)
                sf = SpecFile(filename)
                s1 = sf[0]
                fwhm = s1.scan_header_dict["UFWHM"]
                FWHM_0.append(float(fwhm))
            except:
                FWHM_0.append(0)

        FWHM_1 = []
        for index in Index:
            try:
                filename = "sizes_slit%g_%s/case7keV_%s_spectral_density_%03d.dat" % (APERTURE[1]*1e6, direction, direction, index)
                # filename = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/sizes_slit40.3_h/case7keV_h_spectral_density_000.dat"
                print(">>>>>>>>> opening file: ",filename)
                sf = SpecFile(filename)
                s1 = sf[0]
                fwhm = s1.scan_header_dict["UFWHM"]
                FWHM_1.append(float(fwhm))
            except:
                FWHM_1.append(0)

        FWHM_2 = []
        for index in Index:
            try:
                filename = "sizes_slit%g_%s/case7keV_%s_spectral_density_%03d.dat" % (APERTURE[2]*1e6, direction, direction, index)
                # filename = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/sizes_slit40.3_h/case7keV_h_spectral_density_000.dat"
                print(">>>>>>>>> opening file: ",filename)
                sf = SpecFile(filename)
                s1 = sf[0]
                fwhm = s1.scan_header_dict["UFWHM"]
                FWHM_2.append(float(fwhm))
            except:
                FWHM_2.append(0)

        FWHM_3 = []
        for index in Index:
            try:
                filename = "sizes_slit%g_%s/case7keV_%s_spectral_density_%03d.dat" % (APERTURE[3]*1e6, direction, direction, index)
                # filename = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/sizes_slit40.3_h/case7keV_h_spectral_density_000.dat"
                print(">>>>>>>>> opening file: ",filename)
                sf = SpecFile(filename)
                s1 = sf[0]
                fwhm = s1.scan_header_dict["UFWHM"]
                FWHM_3.append(float(fwhm))
            except:
                FWHM_3.append(0)



        F2theory1 = []
        F2theory2 = []
        Msource_at_id = []
        Msource_at_slit = []

        for index in Index:
            ff_source_at_id, mm_source_at_id = get_f2(f1=F1_0[index],
                                     position_source=0.0,  # source at source
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=False)

            ff_source_at_slit, mm_source_at_slit = get_f2(f1=F1_0[index],
                                     position_source=36.0, #0.0,  # source at slit
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=False)

            Msource_at_id.append(mm_source_at_id)
            Msource_at_slit.append(mm_source_at_slit)

            F2theory1.append(ff_source_at_id)
            F2theory2.append(ff_source_at_slit)


        # print("F:", len(F1), len(F2))


        F2theory1smooth = numpy.array(F2theory1)
        F2theory2smooth = numpy.array(F2theory2)
        Msource_at_id = numpy.array(Msource_at_id)
        Msource_at_slit = numpy.array(Msource_at_slit)


        # fig, ax = plot(
        #     numpy.array(F1_0), numpy.array(F2_0),
        #     numpy.array(F1_1), numpy.array(F2_1),
        #     numpy.array(F1_2), numpy.array(F2_2),
        #     numpy.array(F1_3), numpy.array(F2_3),
        #     numpy.array(F1_0), F2theory1smooth,
        #     numpy.array(F1_0), F2theory2smooth,
        #     legend=["Slit %g um" % (1e6 * APERTURE[0]),
        #             "Slit %g um" % (1e6 * APERTURE[1]),
        #             "Slit %g um" % (1e6 * APERTURE[1]),
        #             "Slit %g um" % (1e6 * APERTURE[2]),
        #             "Geometrical Optics smoothed (source at ID)","Geometrical Optics (source at slit)"],
        #     marker=[None,None,None,None,None,None],
        #     color=['r', 'b', 'g', 'c', 'k', 'k'],
        #     linestyle=[None, None, None, None, '--', ':'],
        #     yrange=[0,60],
        #     xtitle="F1 [m]", ytitle="F2 [m]", title="%s trajectories" % direction,
        #     show=0)
        #
        # ax.xaxis.grid()
        # ax.yaxis.grid()


        fig2, ax2 = plot(
                numpy.array(F1_0), numpy.array(FWHM_0),
                numpy.array(F1_1), numpy.array(FWHM_1),
                numpy.array(F1_2), numpy.array(FWHM_2),
                numpy.array(F1_3), numpy.array(FWHM_3),
                numpy.array(F1_0), Msource_at_id * sourcesize * 1e6,
                numpy.array(F1_0), Msource_at_slit * APERTURE[0] * 1e6,
                marker=[None,None,None,None,None,None],
                color=['r', 'b', 'g', 'c', 'k', 'k'],
                linestyle=[None, None, None, None, '--', ':'],
                yrange=yrange, # numpy.array(FWHM_0).max()*1.1],
                legend=["Slit %g um" % (1e6 * APERTURE[0]),\
                "Slit %g um" % (1e6 * APERTURE[1]), \
                "Slit %g um" % (1e6 * APERTURE[2]), \
                "Slit %g um" % (1e6 * APERTURE[3]), \
                "Analytical (source at ID)","Analytical (source at slit)"],
                xtitle="F1 [m]", ytitle="FWHM [um]", title="%s Sizes" % (direction),
                ylog=1,
                show=0)

        ax2.xaxis.grid()
        ax2.yaxis.grid()

        if graphic_filename is not None:
            plt.savefig(graphic_filename)
            print("File written to disk: %s" % graphic_filename)

        plt.show()

if __name__ == "__main__":


    # plot_one_case('h',0)
    plot_all_cases('h', graphic_filename='sizes_h.eps')
    plot_all_cases('v', graphic_filename='sizes_v.eps')

