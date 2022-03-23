
import numpy
from srxraylib.plot.gol import plot_image, plot_image_with_histograms, plot, plot_show
from tools import fwhm, load_SRW_SD

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
    # SRW
    #
    SD, x, y = plotSD("/nobackup/gurb1/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/intensity7.0keV_15k_ME_intensity.dat",
                      do_plot=0)

    #
    # wofry1d
    #
    h = numpy.loadtxt("spectral_density_data/case1h_wofry_spectral_density.dat")
    v = numpy.loadtxt("spectral_density_data/case1v_wofry_spectral_density.dat")

    # plot(v[:,0], v[:,1],
    #      h[:,0], h[:,1],
    #      title="", legend=[
    #         "V profile (fwhm: %d) " % (fwhm(v[:,0], v[:,1],)),
    #         "H profile (fwhm: %d) " % (fwhm(h[:,0], h[:,1],)),
    #     ],
    #      )



    #
    # comsyl
    #

    from comsyl_propagate36m import add_modes

    comsyl_i, comsyl_x, comsyl_y = add_modes()

    # print(">>>>>>>>>>>>>>>", comsyl_i.shape, comsyl_x.shape, comsyl_y.shape)
    # plot_image(comsyl_i, comsyl_x * 1e6, comsyl_y * 1e6,
    #            title="[COMSYL]",
    #            xtitle="",
    #            ytitle="",
    #            show=1,
    #            )
    #
    # Hybrid
    #


    #
    # wofry1d
    #
    hybrid_h = numpy.loadtxt("spectral_density_data/hybrid36_h.txt")
    hybrid_v = numpy.loadtxt("spectral_density_data/hybrid36_v.txt")

    # plot(hybrid_v[:,0], hybrid_v[:,1],
    #      hybrid_h[:,0], hybrid_h[:,1],
    #      title="", legend=[
    #         "hybrid V profile (fwhm: %d) " % (fwhm(hybrid_v[:,0], hybrid_v[:,1],)),
    #         "hybrid H profile (fwhm: %d) " % (fwhm(hybrid_h[:,0], hybrid_h[:,1],)),
    #     ],
    #      )

    #
    # all
    #

    plot(
         v[:,0], v[:,1] / (v[:,1] ).max(),
         h[:,0], h[:,1] / (h[:,1] ).max(),
         1e6 * y, SD[SD.shape[0] // 2, :] / (SD[SD.shape[0] // 2, :]).max(),
         1e6 * x, SD[:, SD.shape[1] // 2] / (SD[:, SD.shape[1] // 2]).max(),
         1e6 * comsyl_y, comsyl_i[comsyl_i.shape[0] // 2, :] / (comsyl_i[comsyl_i.shape[0] // 2, :]).max(),
         1e6 * comsyl_x, comsyl_i[:, comsyl_i.shape[1] // 2] / (comsyl_i[:, comsyl_i.shape[1] // 2]).max(),
         hybrid_v[:, 0], hybrid_v[:, 1] / (hybrid_v[:, 1]).max(),
         hybrid_h[:, 0], hybrid_h[:, 1] / (hybrid_h[:, 1]).max(),
         title="", legend=[
            "WOFRY1D V profile (fwhm: %d) " % (fwhm(v[:,0], v[:,1],)),
            "WOFRY1D H profile (fwhm: %d) " % (fwhm(h[:,0], h[:,1],)),
            "SRW V profile (fwhm: %d) " % (fwhm(1e6 * y, SD[SD.shape[0] // 2, :])),
            "SRW H profile (fwhm: %d) " % (fwhm(1e6 * x, SD[:, SD.shape[1] // 2])),
            "COMSYL V profile (fwhm: %d) " % (fwhm(1e6 * comsyl_y, comsyl_i[comsyl_i.shape[0] // 2, :] / (comsyl_i[comsyl_i.shape[0] // 2, :]).max(),)),
            "COMSYL H profile (fwhm: %d) " % (fwhm(1e6 * comsyl_x, comsyl_i[:, comsyl_i.shape[1] // 2] / (comsyl_i[:, comsyl_i.shape[1] // 2]).max(),)),
            "HYBRID V profile (fwhm: %d) " % (fwhm(hybrid_v[:, 0], hybrid_v[:, 1], )),
            "HYBRID H profile (fwhm: %d) " % (fwhm(hybrid_h[:, 0], hybrid_h[:, 1], )),
        ],
         color=['red','red','blue','blue','green','green','k','k'],
         linestyle=[None,":",None,":",None,":",None,":"],
         )
