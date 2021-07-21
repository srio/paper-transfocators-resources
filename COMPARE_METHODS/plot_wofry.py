import numpy
import h5py
from srxraylib.plot.gol import plot_image_with_histograms, plot_show
import matplotlib.pylab as plt
from barc4plots.barc4plots import Image2Plot, ESRF_colors_2D
from oasys.util.oasys_util import get_fwhm
import xraylib
from silx.io.specfile import SpecFile
from srxraylib.plot.gol import plot, plot_image


if __name__ == "__main__":
    fontsize = 20
    # fontsize_legend = 22
    # matplotlib.rc('xtick', labelsize=fontsize)
    # matplotlib.rc('ytick', labelsize=fontsize)
    # params = {'legend.fontsize':     fontsize,
    #           'legend.handlelength': fontsize // 20}
    # plt.rcParams.update(params)


    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'figure.autolayout': True})

    plt.rc('font', size=fontsize)
    # plt.rc('axes', titlesize=fontsize)
    plt.rc('axes', labelsize=fontsize * 8 // 10)
    plt.rc('xtick', labelsize=fontsize * 8 // 10)
    plt.rc('ytick', labelsize=fontsize * 8 // 10)
    # plt.tight_layout()





# h/v     slit    F1        F2        R1         R2       FWHM     idx
# h &      40.3 & 46.1 &     26.5 &     641.9 &     369.5 &     8.6 &     86 \\
# v &      227.0 & 85.2 &     27.8 &     1187.4 &     387.3 &     8.8 &     168 \\

# h &      40.3 & 25.1 &     21.3 &     349.1 &     296.3 &     42.0 &     42 \\
# v &      227.0 & 42.2 &     55.6 &     588.6 &     775.3 &     32.3 &     78 \\

# """
# h/v     slit    F1        F2        R1         R2       FWHM     idx
# h &      85.1 & 46.1 &     31.8 &     641.9 &     443.7 &     35.5 &     86 \\
# v &      506.7 & 85.2 &     27.8 &     1187.4 &     387.6 &     5.3 &     168 \\

# h &      85.1 & 25.1 &     20.7 &     349.1 &     288.7 &     24.0 &     42 \\
# v &      506.7 & 42.2 &     55.7 &     588.6 &     776.0 &     134.3 &     78 \\
# """

    print("h/v     slit    F1        F2        R1         R2       FWHM     idx")

    A0_H = []  # abscissas
    A1_H = []  # profile H
    A1N_H = []  # normalized profile H
    A0_V = []  # abscissas
    A1_V = []  # profile H
    A1N_V = []  # normalized profile H


    # LEGEND = []

    Cases       = [1,2,3,4]
    Apertures_h = [40.3e-6, 40.3e-6, 85.1e-6, 85.1e-6]
    Apertures_v = [227.0e-6, 227.0e-6, 506.7e-6, 506.7e-6]
    Selected_f1_h = [46.1, 25.1, 46.1, 25.1]
    Selected_f1_v = [15.0, 42.2, 85.2, 42.2]
    RangeH = [[-50,50],[-50,50],[-50,50],[-50,50]]
    RangeV = [[-200, 200], [-200, 200], [-200, 200], [-200, 200]]


    # for direction in ['h','v']:
    #     # direction = 'v'
    #
    #
    #     if direction == 'h':
    #         APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]
    #         Selected_f1 = [46, 25]
    #
    #         aperture = APERTURE[0]
    #         Cases = [1,2]
    #
    #         # aperture = APERTURE[1]
    #         # Cases = [3,4]
    #     else:
    #         APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]
    #         Selected_f1 = [85, 42]
    #
    #         aperture = APERTURE[1]
    #         Cases = [1,2]
    #
    #         # aperture = APERTURE[2]
    #         # Cases = [3,4]

    for i in range(4):
        for direction in ['h', 'v']:

            print(">>>>>>>>>>>>>>>>>>>", i, direction)
            fileroot = "case7keV_%s" % direction
            if direction == 'h':
                aperture = Apertures_h[i]
            else:
                aperture = Apertures_v[i]

            subdirectory = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV//sizes_slit%g_%s" % (1e6 * aperture, direction)
            a = numpy.loadtxt("/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (1e6 * aperture, direction))
            F1 = a[:, 0].copy()
            F2 = a[:, 1].copy()

            R1 = []
            R2 = []
            for ii in range(F1.size):
                xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
                R1.append(F1[ii] * (2 * xrl_delta))
                R2.append(F2[ii] * (2 * xrl_delta))


            # plot(F1, F2, title="case %d; %s  slit: %g um" % (Cases[i], direction, aperture*1e6))




        # for selected_f1 in Selected_f1:
            if direction == 'h':
                selected_f1 = Selected_f1_h[i]
            else:
                selected_f1 = Selected_f1_v[i]

            index = numpy.argmin(  numpy.abs(F1-selected_f1) )


            filename = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/sizes_slit%g_%s/case7keV_%s_spectral_density_%03d.dat" % (
                aperture * 1e6, direction, direction, index)
            sf = SpecFile(filename)
            s1 = sf[0]
            fwhm = float(s1.scan_header_dict["UFWHM"])

            if direction == 'h':
                A0_H.append(s1.data[0,:])
                A1_H.append(s1.data[1,:])
                A1N_H.append(s1.data[1, :] / s1.data[1, :].max())
            else:
                A0_V.append(s1.data[0,:])
                A1_V.append(s1.data[1,:])
                A1N_V.append(s1.data[1, :] / s1.data[1, :].max())


            # LEGEND.append("%s a=%g um" % (direction, 1e6 * aperture))

            print("%s &      %4.1f & %3.1f &     %3.1f &     %3.1f &     %3.1f &     %3.1f &     %d \\\\" % (direction, 1e6*aperture, F1[index], F2[index],
                                                                                      1e6*R1[index], 1e6*R2[index], fwhm, index))


    print(len(A0_H))
    for i in range(4):
        # plot(A0_H[i], A1N_H[i],
        #      A0_V[i], A1N_V[i], title="Case %d" % i)


        img = numpy.outer(A1N_H[i], A1N_V[i])
        x = A0_H[i]
        y = A0_V[i]

        nx, ny = img.shape

        fwhmH, quote, coordinates = get_fwhm(img[:, ny // 2], x)
        fwhmV, quote, coordinates = get_fwhm(img[nx // 2, :], y)


        plot_image_with_histograms(img, x, y,
                        use_profiles_instead_histograms=1,
                        xtitle="X [$\mu$m]", ytitle="Y [$\mu$m]",
                        title="%4.1f x %4.1f $\mu$m$^2$" % (fwhmH, fwhmV),
                        aspect='auto',
                        add_colorbar=0,
                        figsize=(10,10),
                        cmap=ESRF_colors_2D(8),
                        xrange=RangeH[i], yrange=RangeV[i],
                        show=0)

        filepng = "case%d_wofry.png" % Cases[i]
        plt.savefig(filepng)
        print("File written to disk: %s" % filepng)

        plot_show()


        # fig, ax = plot(A0[0], A1N[0],
        #      A0[2], A1N[2],
        #      legend=["case %d h" % Cases[0], "case %d v" % Cases[0]], figsize=[10,10], show=0)
        # ax.xaxis.grid()
        # ax.yaxis.grid()
        # ax.set_xlabel("x [$\mu$ m]", fontsize=20)
        # ax.set_ylabel("Intensity [a.u.]", fontsize=20)
        #
        # filename = "case_%d_profiles.eps" % Cases[0]
        # plt.savefig(filename)
        # print("File written to disk: %s" % filename)
        # plt.show()
    #
    #
    # fig, ax = plot(A0[1], A1N[1],
    #      A0[3], A1N[3],
    #      legend=["case %d h" % Cases[1], "case %d v" % Cases[1]], figsize=[10,10], show=0)
    # ax.xaxis.grid()
    # ax.yaxis.grid()
    # ax.set_xlabel("y [$\mu$ m]", fontsize=20)
    # ax.set_ylabel("Intensity [a.u.]", fontsize=20)
    #
    # filename = "case_%d_profiles.eps" % Cases[1]
    # plt.savefig(filename)
    # print("File written to disk: %s" % filename)
    # plt.show()
    #
    # title="" #case %d" % Cases[0]
    # fig, ax = plot_image( numpy.outer(A1N[0], A1N[2]), A0[0], A0[2], title=title, figsize=[10,10], aspect='auto', show=0)
    # ax.set_xlabel("x [$\mu$ m]", fontsize=20)
    # ax.set_ylabel("y [$\mu$ m]", fontsize=20)
    # filename = "case_%d_image.eps" % Cases[0]
    # plt.savefig(filename)
    # print("File written to disk: %s" % filename)
    # plt.show()
    #
    # title="" #case %d" % Cases[1]
    # fig, ax = plot_image( numpy.outer(A1N[1], A1N[3]), A0[1], A0[3], title=title, figsize=[10,10], aspect='auto', show=0)
    # ax.set_xlabel("x [$\mu$ m]", fontsize=20)
    # ax.set_ylabel("y [$\mu$ m]", fontsize=20)
    # filename = "case_%d_image.eps" % Cases[1]
    # plt.savefig(filename)
    # print("File written to disk: %s" % filename)
    # plt.show()
    #







    # from matplotlib import rcParams
    #
    #
    # dir = "/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/results_configs/"
    #
    # for i in range(1,5):
    #     file = "id18_c0%d_7.0keV_25k_ME_intensity.dat" % i
    #     print("file: ", dir+file)
    #
    #     x, y, img = native_util.load_intensity_file(dir+file)
    #
    #     x *= 1e6
    #     y *= 1e6
    #
    #     (fwhmH, fwhmV) = (0, 0)
    #     nx, ny = img.shape
    #
    #     fwhmH, quote, coordinates = get_fwhm(img[:, ny // 2], x)
    #     fwhmV, quote, coordinates = get_fwhm(img[nx // 2, :], y)
    #
    #     fig, ai, ax, ay = plot_image_with_histograms(img.T,x,y, use_profiles_instead_histograms=1,
    #                                                  xtitle="X [$\mu$m]", ytitle="Y [$\mu$m]",
    #                                                  title="%4.1f x %4.1f $\mu$m$^2$" % (fwhmH, fwhmV),
    #                                                  aspect_ratio='auto',
    #                                                  add_colorbar=0,
    #                                                  figsize=(10,10),
    #                                                  cmap=ESRF_colors_2D(8),
    #                                                  show=0)
    #
    #     filepng = "case%d_srw.png" % i
    #     plt.savefig(filepng)
    #     print("File written to disk: %s" % filepng)
    #
    #     plot_show()


