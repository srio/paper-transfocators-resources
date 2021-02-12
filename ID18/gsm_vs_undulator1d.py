
#
# Import section
#
import numpy

def W(x1,x2):
    delta_x = x2 - x1
    mu = numpy.exp(-delta_x ** 2 / 2 / sigma_xi ** 2)
    s1 = numpy.exp(-x1 ** 2 / 2 / sigma_x ** 2)
    s2 = numpy.exp(-x2 ** 2 / 2 / sigma_x ** 2)
    return mu * numpy.sqrt(s1) * numpy.sqrt(s2)

def get_q(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    return q

def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF

if __name__ == "__main__":

    from srxraylib.plot.gol import plot, plot_image, plot_table

    # electrons
    sigmaxx = 3.01836e-05
    sigmaxpxp = 4.36821e-06

    # GSM
    beta = 0.0922395 #1.151 #0.02 #
    sigma_x = 3.03783e-05
    sigma_xi = beta * sigma_x



    x1 = numpy.linspace(-0.00012, 0.00012, 400)

    X1 = numpy.outer(x1, numpy.ones_like(x1))
    X2 = numpy.outer(numpy.ones_like(x1), x1 )

    CSD_GSM = W(X1,X2)
    indices = numpy.arange(x1.size)

    # plot(x1, W(x1,x1),
    #      x1, CSD_GSM[indices, indices],
    #      x1, numpy.exp(-x1**2/2/sigma_x**2),
    #      title="Spectral density", legend=["SD function", "SD array", "Gaussian with sigma_x"])
    #
    # plot_image(CSD_GSM, x1, x1)

    from wofry.propagator.util.gaussian_schell_model import GaussianSchellModel1D, GaussianSchellModel2D

    nmodes = 100
    gsm = GaussianSchellModel1D(1.0, sigma_x, sigma_xi)

    # plot(x1, y2) # , legend=["numeric","theoretical"] )

    eigenvalues_GSM = numpy.zeros(nmodes)
    eigenvectors_GSM = numpy.zeros((nmodes,x1.size))
    q = get_q(beta)
    lambda0 = numpy.sqrt(numpy.pi) / (gsm.a() + gsm.b() + gsm.c())
    for i in range(nmodes):
        eigenvalues_GSM[i] = lambda0 * q**i
        eigenvectors_GSM[i,:] = gsm.phi(i, x1)

    occupation_GSM = eigenvalues_GSM / eigenvalues_GSM.sum()
    cumulated_occupation_GSM = numpy.cumsum(occupation_GSM)

    indices = numpy.arange(x1.size)
    spectral_density_GSM = CSD_GSM[indices,indices]

    print("CF: ", occupation_GSM[0])
    print(">>>>", eigenvectors_GSM.shape)

    plot(numpy.arange(nmodes), cumulated_occupation_GSM, title="Cumulated occupation")
    plot(x1, spectral_density_GSM, title="Spectral density")
    plot_table(x1, eigenvectors_GSM, title="Eigenvectors")



    #
    # undulator
    #
    from undulator_coherent_mode_decomposition_1d import UndulatorCoherentModeDecomposition1D

    co = UndulatorCoherentModeDecomposition1D(
                                    electron_energy=6.0,
                                    electron_current=0.2,
                                    undulator_period=0.02,
                                    undulator_nperiods=100,
                                    K=1.191085,
                                    photon_energy=10000.,
                                    nsigma=16,
                                    number_of_points=200,
                                    distance_to_screen=100.0,
                                    scan_direction="V",
                                    sigmaxx=sigmaxx,
                                    sigmaxpxp=sigmaxpxp)
    # make calculation
    out = co.calculate()

    CSD_U = out["CSD"]
    abscissas = out["abscissas"]
    eigenvalues_U = out["eigenvalues"]
    eigenvectors_U = out["eigenvectors"]

    occupation_U = eigenvalues_U[0:nmodes] / eigenvalues_U.sum()
    cumulated_occupation_U = numpy.cumsum(occupation_U)
    indices = numpy.arange(abscissas.size)
    spectral_density_U = CSD_U[indices,indices]

    plot(numpy.arange(nmodes), cumulated_occupation_U, title="Cumulated occupation U")
    # plot(abscissas, spectral_density_U, title="Spectral density U")
    plot_table(abscissas, eigenvectors_U[0:10,:], title="Eigenvectors U")

    # restore spectral density from modes
    y = numpy.zeros_like(abscissas, dtype=complex)
    for i in range(200):
        y += eigenvalues_U[i] * numpy.conjugate(eigenvectors_U[i, :]) * eigenvectors_U[i, :]

    y = numpy.real(y)
    plot(x1, spectral_density_GSM / spectral_density_GSM.max(),
         abscissas, spectral_density_U / spectral_density_U.max(),
         abscissas, y / y.max(), legend=["SD_GSM", "SD", "SD From modes"])

    plot(numpy.arange(nmodes), cumulated_occupation_GSM,
         numpy.arange(nmodes), cumulated_occupation_U,
         title="Cumulated occupation", legend=["GSM","U"])

    for i in range(10):
        plot(x1, eigenvectors_GSM[i, :] / (eigenvectors_GSM[i, :]).max(),
             abscissas, eigenvectors_U[i, :] / (eigenvectors_U[i, :]).max(),
             title="Mode: %d" % i)

