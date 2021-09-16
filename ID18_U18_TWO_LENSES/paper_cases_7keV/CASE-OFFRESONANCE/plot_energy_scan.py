import numpy
from srxraylib.plot.gol import plot


Energy = numpy.linspace(6900,7100,100)

cumulated = 0.0
spectrum = numpy.zeros_like(Energy)
cf = numpy.zeros_like(Energy)


plot_from_raw_files = False
what = "source"
what = "sample"

if plot_from_raw_files:
    dir = "case1h/"
    e0 = 7000
    deltahalf = 1
    for i,energy in enumerate(Energy):
        if energy > (e0-20000*deltahalf)  and energy < (e0+10000*deltahalf):
            try:
                sim = numpy.loadtxt("%s/%s_spectral_density_%4d.dat" % (dir,what, energy),
                                    skiprows=4)
                current = sim[:,1]
                cumulated = cumulated + current
                spectrum[i] = current.sum() * (sim[1, 0] - sim[0, 0])

                if False:
                    plot(
                        sim[:, 0], current * 10,
                        sim[:, 0], cumulated,
                        xrange=[-1200, 1200],
                        legend=["monochromatic (E=%d eV) x 10" % energy, "pink"],
                        color=["r", "pink"],
                        show=1, ylog=1,
                    )

                cfi = numpy.loadtxt("%s/%s_occupation_%4d.dat" % (dir, what, energy),
                                    skiprows=3)
                cf[i] = cfi[0,1]
            except:
                break


    e00 = 7001
    sim = numpy.loadtxt("%s/%s_spectral_density_%4d.dat" % (dir, what, e00), skiprows=4)
    current = sim[:,1]

    plot(
        sim[:, 0], current / current.max(),
        sim[:, 0], cumulated/ cumulated.max(),
        xrange=[-40,40],
        legend=["monochromatic (E=%d eV) x 10" % e00, "pink"],
        color=["r","pink"],
        show=0,
        )

    plot(Energy, spectrum, show=0)
    plot(Energy, cf, title="CF")


    filename = "case1h_intensity_at_%s.dat" % what
    f = open(filename, "w")
    for i in range(cumulated.size):
        f.write("%g  %g  %g\n" % (sim[i, 0], cumulated[i], current[i]))
    f.close()
    print("File %s written to disk" % filename)

    filename = "case1h_energy_scan_at_%s.dat" % what
    f = open(filename, "w")
    for i in range(Energy.size):
        f.write("%g  %g  %g\n" % (Energy[i], spectrum[i], cf[i]))
    f.close()
    print("File %s written to disk" % filename)

else:
    #
    # replot from colected files
    #

    sim = numpy.loadtxt("case1h_intensity_at_%s.dat" % what)
    e00 = 7001
    plot(
        sim[:, 0], sim[:, 2] / sim[:, 2].max(),
        sim[:, 0], sim[:, 1] / sim[:, 1].max(),
        xrange=[-40,40],
        legend=["monochromatic (E=%d eV)" % e00, "pink"],
        color=["r","pink"],
        show=0,
        )


    ene = numpy.loadtxt("case1h_energy_scan_at_%s.dat" % what)
    plot(ene[:,0], ene[:,1], show=0)
    plot(ene[:,0], ene[:,2], title="CF")

