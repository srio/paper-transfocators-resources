import numpy
from srxraylib.plot.gol import plot_image_with_histograms, plot_show
from oasys.util.oasys_util import get_fwhm
import Shadow
from barc4plots.barc4plots import Image2Plot, ESRF_colors_2D
import matplotlib.pylab as plt
# from srxraylib.plot.gol import set_qt
# set_qt()

def get_tkt_fwhm(tkt):

    x = tkt["bin_h_center"]
    y = tkt["bin_v_center"]
    hx = tkt["histogram_h"]
    hy = tkt["histogram_v"]

    fwhm_x, _, _ = get_fwhm(hx, x)
    fwhm_y, _, _ = get_fwhm(hy, y)
    return fwhm_x, fwhm_y


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

    dir = "/users/srio/Oasys/to_Manolo/shadow_images_star_files"

    for i in range(1,5):

        # filename = "%s/%dh_%dv/star.08" % (dir, i, i)
        filename = "%s/%dh_%dv/screen.0801" % (dir, i, i)
        beam = Shadow.Beam()
        beam.load(filename)

        w = 300e-6
        tkt = beam.histo2(1, 3, ref=23, nolost=1, nbins=401, xrange=[-w / 2, w / 2], yrange=[-w / 2, w / 2])

        hfactor = 1e6
        vfactor = 1e6

        img = tkt["histogram"]
        nx, ny = img.shape

        x = numpy.linspace(hfactor * tkt["xrange"][0], hfactor * tkt["xrange"][1], nx)
        y = numpy.linspace(vfactor * tkt["yrange"][0], vfactor * tkt["yrange"][1], ny)


        # ah = tkt["bin_h_edges"]
        # ah = tkt["bin_v_edges"]


        fwhmH, quote, coordinates = get_fwhm(img[:, (ny // 2)], x)
        fwhmV, quote, coordinates = get_fwhm(img[(nx // 2), :], y)


        # width = 1
        #
        # fwhmH, quote, coordinates = get_fwhm(img[:, (ny // 2 - width):(ny // 2 + width)], x)
        # fwhmV, quote, coordinates = get_fwhm(img[(nx // 2 - width):(nx // 2 + width), :], y)


        plot_image_with_histograms(img, x, y,
                        use_profiles_instead_histograms=0,
                        xtitle="X [$\mu$m]", ytitle="Y [$\mu$m]",
                        title="%4.1f x %4.1f $\mu$m$^2$" % (fwhmH, fwhmV),
                        aspect='auto',
                        add_colorbar=0,
                        figsize=(10,10),
                        cmap=ESRF_colors_2D(8),
                        xrange=RangeH[i-1], yrange=RangeV[i-1],
                        show=0)

        filepng = "case%d_hybrid.png" % (i)
        plt.savefig(filepng)
        print("File written to disk: %s" % filepng)

        plot_show()