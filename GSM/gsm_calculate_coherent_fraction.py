import numpy

def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF


if __name__ == "__main__":
    from srxraylib.plot.gol import plot
    #
    # coherent fraction
    #

    betas = numpy.linspace(0.,10.0, 100)

    CF = get_coherent_fraction_exact(betas)
    f = plot(betas,CF,xtitle=r'$\mathrm{\beta}$',ytitle="CF",show=True,figsize=(10,7.5))


