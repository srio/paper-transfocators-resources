
import numpy
from srxraylib.plot.gol import plot_image, plot_image_with_histograms, plot, plot_show

from tools import fwhm, load_SRW_SD, plot_four_images_with_histograms

import pylab as plt



def plotSD(intensity_file_name, do_plot=0):


    SD, x, y = load_SRW_SD(intensity_file_name)




    # plot_image_with_histograms(intensity, 1e6 * x, 1e6 * y,
    #                            title="[SRW]",
    #                            xtitle="",
    #                            ytitle="",
    #                            show=0,
    #                            use_profiles_instead_histograms=True,
    #                            xrange=range_limits,
    #                            yrange=range_limits,
    #                            )

    if do_plot:
        plot_image(SD, x * 1e6, y * 1e6,
                            title="[SRW]",
                            xtitle="",
                            ytitle="",
                            show=0,
                            # xrange=range_limits,
                            # yrange=range_limits,
                   )


        plot(
             1e6*y, SD[SD.shape[0] // 2, :],
             1e6*x, SD[:, SD.shape[0] // 2],
             title="", legend=[
                "V profile (fwhm: %d) " %    (fwhm(1e6*y, SD[SD.shape[0] // 2, :])),
                "H profile (fwhm: %d) " %    (fwhm(1e6*x, SD[:, SD.shape[0] // 2])),
            ],
             xtitle="",
             ytitle="",
             show=0)


        plot_show()

    return SD, x, y

if __name__ == "__main__":
    INT = []
    X = []
    Y = []
    TITLE = []




    #
    # wofry1d
    #
    h = numpy.loadtxt("spectral_density_data/case1h_wofry_spectral_density.dat")
    v = numpy.loadtxt("spectral_density_data/case1v_wofry_spectral_density.dat")

    hx = h[:,0]
    hy = v[:,0]
    hxi = h[:,1]
    hyi = v[:,1]
    # plot(hy,hyi,
    #      hx, hxi,
    #      title="", legend=[
    #         "V profile (fwhm: %d) " % (fwhm(hy, hyi)),
    #         "H profile (fwhm: %d) " % (fwhm(hx, hxi)),
    #     ],
    #      )

    INT.append( numpy.outer( hxi, hyi))
    X.append(hx)
    Y.append(hy)
    TITLE.append("WOFRY1D")


    #
    # comsyl
    #

    from comsyl_propagate36m import add_modes

    comsyl_i, comsyl_x, comsyl_y = add_modes()

    INT.append(comsyl_i)
    X.append(comsyl_x*1e6)
    Y.append(comsyl_y*1e6)
    TITLE.append("COMSYL")




    #
    # SRW
    #
    SD, x, y = plotSD("/nobackup/gurb1/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/intensity7.0keV_15k_ME_intensity.dat",
                      do_plot=0)


    INT.append(SD)
    X.append(x*1e6)
    Y.append(y*1e6)
    TITLE.append("SRW")



    #
    # hybrid
    #
    import h5py
    filename = "spectral_density_data/hybrid_xy_plot.hdf5"
    f = h5py.File(filename,'r')

    x = f['coordinates/X'][:]
    y = f['coordinates/Y'][:]
    ii = f['xy_plots/last_plot/23: Total Intensity = |Eσ|² + |Eπ|²'][:]
    f.close()

    INT.append(ii)
    X.append(x*1e6)
    Y.append(y*1e6)
    TITLE.append("HYBRID")


    ff = plot_four_images_with_histograms(INT, X, Y,
                               title_list=TITLE,
                               xtitle="x [$\mu$m]",
                               ytitle="y [$\mu$m]",
                               show=0,
                               use_profiles_instead_histograms=True,
                               xrange=[-1000,1000],
                               yrange=[-1000,1000],
                               cmap=plt.cm.hsv,
                               )

    plt.savefig("plot_2D_spectral_density_36m.png")
    plot_show()
    print("File written to disk: plot_2D_spectral_density_36m.png")
