import numpy
from srxraylib.plot.gol import plot, plot_show




if __name__ == "__main__":


    #
    # New Journal of Physics 12 (2010) 035004  I A Vartanyants 1 and A Singer
    #



    z = numpy.linspace(0., 100.0, 100)
    CF = 0.9
    wavelength = 1e-10

    sigma = 10e-6
    beta = CF / numpy.sqrt(1 - CF)
    mu = beta * sigma
    sigmap = (wavelength/4/numpy.pi) / sigma

    # eq 18
    delta2 = 1 / ( 1 / (4*sigma**2) + 1 / mu**2)
    delta = numpy.sqrt(delta2)


    k = 2 * numpy.pi / wavelength

    # eq 22
    zeff = k * sigma * delta

    # eq 20
    Deltaz = numpy.sqrt(1 + (z/zeff)**2)



    if False:
        plot(z,Deltaz,
             z, z/z * numpy.sqrt(2),
             xtitle="z [m]",
             legend=["D(z)","D(zeff)"])

    #
    # calculation of propagated mu
    #



    Mu0 = delta * Deltaz # eq 19
    Mu0 = mu * Deltaz  # eq 19 CORRECTED!!!!!!!!!!!!!!!!!!!!!!


    #
    # Sigmap = 1 / (2*k*mu) * numpy.sqrt(4+beta**2) # eq 25
    theta_Chi = 1 / (2*k*sigma) * numpy.sqrt(4+beta**2)
    Mu1 = mu * Deltaz # eq 28.1
    Mu2 = numpy.sqrt( mu**2 + (theta_Chi*z)**2) # eq 28.2

    Sigma = sigma * Deltaz # eq 24
    Mu3 = beta * Sigma # see paragraph after eq 33


    print('mu=',mu)
    print('delta=', delta)
    print('zeff=', zeff)

    # Mu4 = numpy.sqrt(mu**2 + (sigmap * z)**2)
    plot(z, 1e6 * Mu0,
         z, 1e6 * Mu1,
         z, 1e6 * Mu2,
         z, 1e6 * Mu3,
         xtitle="z [m]",
         legend=["Mu0'","Mu1", "Mu2", "Mu3"])
