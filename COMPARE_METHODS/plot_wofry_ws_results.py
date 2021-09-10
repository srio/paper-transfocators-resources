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



    RangeH = [[-50,50],[-50,50],[-50,50],[-50,50]]
    RangeV = [[-200, 200], [-200, 200], [-200, 200], [-200, 200]]

    dir = "/users/srio/Oasys/"
    dir = "/users/srio/OASYS1.2/paper-transfocators-resources/COMPARE_METHODS/data_wofry/"

    for i in range(1,5):

        filename_h = dir + "case%d%s_wofry_spectral_density.dat" % (i, 'h')
        filename_v = dir + "case%d%s_wofry_spectral_density.dat" % (i, 'v')
        ah = numpy.loadtxt(filename_h,skiprows=3)
        av = numpy.loadtxt(filename_v, skiprows=3)
        print(ah.shape)
        A0_H = ah[:, 0]
        A1N_H = ah[:, 1]
        A0_V = av[:, 0]
        A1N_V = av[:, 1]
        img = numpy.outer(A1N_H, A1N_V)
        x = A0_H
        y = A0_V

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
                        cmap=ESRF_colors_2D(5),
                        xrange=RangeH[i-1], yrange=RangeV[i-1],
                        show=0)

        filepng = "case%d_wofry_ws_results.png" % (i)
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


