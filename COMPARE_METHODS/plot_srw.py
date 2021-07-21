import numpy
import h5py
from srxraylib.plot.gol import plot_image_with_histograms, plot_show
import matplotlib.pylab as plt
from barc4plots.barc4plots import Image2Plot, ESRF_colors_2D
from orangecontrib.srw.widgets.native.util import native_util
from oasys.util.oasys_util import get_fwhm


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

# from matplotlib import rcParams


dir = "/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/results_configs/"

RangeH = [[-50, 50], [-50, 50], [-50, 50], [-50, 50]]
RangeV = [[-200, 200], [-200, 200], [-200, 200], [-200, 200]]

for i in range(1,5):
    file = "id18_c0%d_7.0keV_25k_ME_intensity.dat" % i
    print("file: ", dir+file)

    x, y, img = native_util.load_intensity_file(dir+file)

    x *= 1e6
    y *= 1e6

    nx, ny = img.shape

    fwhmH, quote, coordinates = get_fwhm(img[:, ny // 2], x)
    fwhmV, quote, coordinates = get_fwhm(img[nx // 2, :], y)

    fig, ai, ax, ay = plot_image_with_histograms(img,x,y, use_profiles_instead_histograms=1,
                                                 xtitle="X [$\mu$m]", ytitle="Y [$\mu$m]",
                                                 title="%4.1f x %4.1f $\mu$m$^2$" % (fwhmH, fwhmV),
                                                 aspect='auto',
                                                 add_colorbar=0,
                                                 figsize=(10,10),
                                                 cmap=ESRF_colors_2D(8),
                                                 xrange=RangeH[i-1], yrange=RangeV[i-1],
                                                 show=0)

    filepng = "case%d_srw.png" % i
    plt.savefig(filepng)
    print("File written to disk: %s" % filepng)

    plot_show()


