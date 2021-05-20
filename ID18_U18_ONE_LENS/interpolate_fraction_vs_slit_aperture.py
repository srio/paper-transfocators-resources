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



    energy_in_keV = 7 # 30 # 15 # 7



    source = "UND"
    aperture = "RECTANGULAR"


    #
    #
    #

    CF_target = 0.90

    # sources = ["UND", "GSM"]
    # apertures = ["RECTANGULAR","GAUSSIAN"]
    # for aperture in apertures:
    #     for source in sources:


    # plot(sl, cfh, sl, cfv)

    # print(">>>>", numpy.all(numpy.diff(cfh) > 0))


    for CF_target in [0.9, 0.8, 0.7]:
        print("\n")
        for energy_in_keV in [7,15,30]:
            subdirectory = "DataCF%d" % energy_in_keV
            filename = "%s/coherent_fraction_vs_slit_source_%s_aperture_%s.dat" % (subdirectory, source, aperture)
            a = numpy.loadtxt(filename)
            sl = a[:, 0]
            cfh = a[:, 1]
            cfv = a[:, 2]

            slh = numpy.interp(CF_target, numpy.flip(cfh), numpy.flip(sl))
            slv = numpy.interp(CF_target, numpy.flip(cfv), numpy.flip(sl))

            print("For energy %2d keV, and CF: %g I need slit %5.1f X %5.1f um2" %(energy_in_keV, CF_target, 1e6*slh, 1e6*slv))
