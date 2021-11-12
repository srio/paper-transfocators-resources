#
# Import section
#
import numpy

#
# MAIN========================
#
if __name__ == "__main__":
    import matplotlib.pylab as plt
    from srxraylib.plot.gol import plot, plot_show, set_qt
    set_qt()

    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'figure.autolayout': True})

    params = {'legend.fontsize': 23,
              'legend.frameon': False,
              'legend.handlelength': 2,
              # 'axes.titlesize' : 24,
              'axes.labelsize': 40,
              'lines.linewidth': 3,
              'lines.markersize': 10,
              'xtick.labelsize': 40,
              'ytick.labelsize': 40,
              'legend.loc': 'upper right',  # 'best'
              # 'grid.color':       'r',
              # 'grid.linestyle':   '-',
              # 'grid.linewidth':     2,
              }
    plt.rcParams.update(params)



    energy_in_keV = 7 # 30 # 15 # 7



    sources = ["UND"] #, "GSM"]
    apertures = ["RECTANGULAR"] #,"GAUSSIAN"]
    subdirectory = "DataCF%d" % energy_in_keV
    #
    #
    #

    #
    # plots
    #
    #
    # paper
    #

    aperture = apertures[0]
    source = sources[0]
    # b7_1 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_GSM_aperture_%s.dat" % ("DataCF7", aperture))
    b7_2 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_UND_aperture_%s.dat" % ("DataCF7", aperture))

    # b15_1 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_GSM_aperture_%s.dat" % ("DataCF15", aperture))
    b15_2 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_UND_aperture_%s.dat" % ("DataCF15", aperture))

    # b30_1 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_GSM_aperture_%s.dat" % ("DataCF30", aperture))
    b30_2 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_UND_aperture_%s.dat" % ("DataCF30", aperture))


    #
    # smooth
    #
    # from scipy.interpolate import make_interp_spline, BSpline
    # T = 1e6 * b2[:, 0]
    # power = b2[:, 1]
    # xnew = numpy.linspace(T.min(), T.max(), T.size)
    # spl = make_interp_spline(T, power, k=1)  # type: BSpline
    # power_smooth = spl(xnew)


    from scipy.signal import savgol_filter
    w7 = savgol_filter(b7_2[:, 1], 5, 2)
    w15 = savgol_filter(b15_2[:, 1], 5, 2)
    w30 = savgol_filter(b30_2[:, 1], 5, 2)



    g = plot(
        1e6 * b7_2[:, 0], w7, #b2[:, 1],
        1e6 * b7_2[:, 0], b7_2[:, 2],
        1e6 * b15_2[:, 0], w15,  # b2[:, 1],
        1e6 * b15_2[:, 0], b15_2[:, 2],
        1e6 * b30_2[:, 0], w30,  # b2[:, 1],
        1e6 * b30_2[:, 0], b30_2[:, 2],
        legend=['7 keV Horizontal', '7 keV Vertical',
                '15 keV Horizontal', '15 keV Vertical',
                '30 keV Horizontal', '30 keV Vertical',
                ],
        figsize=(15,10),
        xtitle=r'Slit aperture [$\mu$m]', ytitle="Coherent Fraction",
        color=['green', 'green', 'blue', 'blue', 'red', 'red'],
        linestyle=[None, '--', None, '--', None, '--',],
        xlog=False, xrange=[0.,1000], yrange=[0, 1.01], show=False)

    g[1].grid()

    plt.yticks(numpy.arange(0, 1.1, step=0.2))

    def aa(x):
        return x / 565
    def invaa(x):
        return 565 * x

    plt.savefig("cf_vs_aperture.eps")

    plot_show()