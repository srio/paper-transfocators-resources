import numpy

from srxraylib.plot.gol import plot, plot_show


if __name__ == "__main__":


    #
    # coherent fraction
    #

    CF = numpy.linspace(0., 1.0, 100)
    sigmas = [1e-6, 10e-6, 100e-6]

    RATIO = []


    for sigma in sigmas:
        a = 1 / 4 / sigma ** 2

        betas = CF / numpy.sqrt(1 - CF)

        mu = betas * sigma

        b = 1 / 2 / mu ** 2

        c = numpy.sqrt(a ** 2 + 2 * a * b)

        sigma_c = numpy.sqrt(1 / 2 / c)

        ratio = sigma_c / mu

        RATIO.append(ratio)

    f, ax = plot(
        CF, numpy.array(RATIO[0]),
        CF, numpy.array(RATIO[1]),
        CF, numpy.array(RATIO[2]),
        xtitle='CF', ytitle="sigma_c/mu", show=False, figsize=(10, 7.5),
        legend=["sigma_x= 1 um","sigma_x= 10 um","sigma_x= 100 um"])

    ax.grid()

    plot_show()



