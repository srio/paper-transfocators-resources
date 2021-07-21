import numpy
import h5py
from srxraylib.plot.gol import plot_image_with_histograms, plot_show
import matplotlib.pylab as plt
from barc4plots.barc4plots import Image2Plot, ESRF_colors_2D

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


dir = "/users/srio/Oasys/"

RangeH = [[-50, 50], [-50, 50], [-50, 50], [-50, 50]]
RangeV = [[-200, 200], [-200, 200], [-200, 200], [-200, 200]]

for i in range(1,5):
    file = "case%d_comsyl.h5" % i
    print("file: ", dir+file)
    f = h5py.File(dir+file, 'r')

    img = f["accumulated_intensity/intensity/intensity_accumulated"][()].T
    x = f["accumulated_intensity/intensity/x"][()]
    y = f["accumulated_intensity/intensity/y"][()]
    fwhmH = f["accumulated_intensity/profileH/fwhm"][()]
    fwhmV = f["accumulated_intensity/profileV/fwhm"][()]

    f.close()


    print(fwhmV)
    fig, ai, ax, ay = plot_image_with_histograms(img,x,y, use_profiles_instead_histograms=1,
                                                 xtitle="X [$\mu$m]", ytitle="Y [$\mu$m]",
                                                 title="%4.1f x %4.1f $\mu$m$^2$" % (fwhmH, fwhmV),
                                                 aspect='auto',
                                                 add_colorbar=0,
                                                 figsize=(10,10),
                                                 cmap=ESRF_colors_2D(8),
                                                 xrange=RangeH[i - 1], yrange=RangeV[i - 1],
                                                 show=0)

    filepng = "case%d_comsyl.png" % i
    plt.savefig(filepng)
    print("File written to disk: %s" % filepng)


    plot_show()


