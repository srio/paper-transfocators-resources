import numpy
from srxraylib.plot.gol import plot, plot_show




if __name__ == "__main__":


    #
    # coherent fraction
    #

    CF = numpy.linspace(0., 1.0, 100)
    sigmas = [1e-6, 10e-6, 100e-6]


    sigma = sigmas[0]



    betas = CF / numpy.sqrt(1 - CF)

    mu = betas * sigma

    # Eq 36 in Santarsiero et al. 106 J. Opt. Soc. Am. A / Vol. 16, No. 1 / January 1999
    M2 = numpy.sqrt(1 + 4 * sigma**2 / mu**2)

    f, ax = plot(
        1/CF, M2, xtitle="1/CF", ytitle="M2",
        )

    f, ax = plot(
        M2, M2 - (2 - CF)/CF , xtitle="M2 - (2 - CF)/CF", ytitle="",
        )





