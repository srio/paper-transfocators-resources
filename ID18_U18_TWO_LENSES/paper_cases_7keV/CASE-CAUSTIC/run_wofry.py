import numpy


#
# MAIN FUNCTION========================
#


def main1(distance=30.0, nmodes=10, directory="xxx", do_plot=True):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(nmodes):
        output_wavefront = run_source(my_mode_index=my_mode_index)
        output_wavefront = run_beamline(output_wavefront, distance=distance)
        tally.append(output_wavefront)

    if do_plot:
        tally.plot_cross_spectral_density(show=1, filename="")
        tally.plot_spectral_density(show=1, filename="")
        tally.plot_occupation(show=1, filename="")

    tally.save_spectral_density(filename="%s/%s_wofry_spectral_density_%g.dat" % (directory, directory, distance))
    # tally.save_occupation(filename="case4v_wofry_occupation.dat")


#
# MAIN========================
#

# from case2h_wofry import run_source, run_beamline
#
# distances = numpy.linspace(10.0, 50.0, 81)
# for distance in distances:
#     main1(distance=distance, nmodes=10, directory="case2h", do_plot=0)
#
#
# from case2v_wofry import run_source, run_beamline
# distances = numpy.linspace(10.0, 50.0, 81)
# for distance in distances:
#     main1(distance=distance, nmodes=10, directory="case2v", do_plot=0)
#
#
#
# from case3h_wofry import run_source, run_beamline
# distances = numpy.linspace(10.0, 50.0, 81)
# for distance in distances:
#     main1(distance=distance, nmodes=10, directory="case3h", do_plot=0)
#
#
# from case3v_wofry import run_source, run_beamline
# distances = numpy.linspace(10.0, 50.0, 81)
# for distance in distances:
#     main1(distance=distance, nmodes=10, directory="case3v", do_plot=0)

# from case4h_wofry import run_source, run_beamline
# distances = numpy.linspace(10.0, 50.0, 81)
# for distance in distances:
#     main1(distance=distance, nmodes=10, directory="case4h", do_plot=0)


# from case2h_wofry import run_source, run_beamline
# distances = numpy.linspace(10.0, 50.0, 81)
# for distance in distances:
#     main1(distance=distance, nmodes=1, directory="case2h0", do_plot=0)

from case4vRefined_wofry import run_source, run_beamline
distances = numpy.linspace(10.0, 50.0, 81)
for distance in distances:
    main1(distance=distance, nmodes=10, directory="case4vRefined", do_plot=0)