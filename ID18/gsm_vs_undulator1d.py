
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

    do_plot = False
    scan_direction = "H"

    if scan_direction == "H":
        # HORIZONTAL ----------------------
        # electrons
        sigmaxx = 3.01836e-05
        sigmaxpxp = 4.36821e-06
        # GSM
        beta = 0.0922395 #0.02 #
        sigma_x = 3.03783e-05
        x1 = numpy.linspace(-0.00012, 0.00012, 400)
        nmodes = 100
        nsigma = 16
    else:
        # VERTICAL ----------------------
        # electrons
        sigmaxx = 3.63641e-06
        sigmaxpxp = 1.37498e-06
        # GSM
        beta = 1.151 #0.02 #
        sigma_x = 5.001e-6
        x1 = numpy.linspace(-0.000050, 0.000050, 400)
        nmodes = 10
        nsigma = 6

    #
    #
    #

    sigma_xi = beta * sigma_x

    X1 = numpy.outer(x1, numpy.ones_like(x1))
    X2 = numpy.outer(numpy.ones_like(x1), x1 )

    CSD_GSM = W(X1,X2)
    indices = numpy.arange(x1.size)

    if do_plot:
        plot(x1, W(x1,x1),
             x1, CSD_GSM[indices, indices],
             x1, numpy.exp(-x1**2/2/sigma_x**2),
             title="Spectral density", legend=["SD function", "SD array", "Gaussian with sigma_x"])


    from wofry.propagator.util.gaussian_schell_model import GaussianSchellModel1D, GaussianSchellModel2D


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

    if do_plot:
        plot(numpy.arange(nmodes), cumulated_occupation_GSM, title="Cumulated occupation GSM")
        plot(x1, spectral_density_GSM, title="Spectral density GSM")
        plot_table(x1, eigenvectors_GSM, title="Eigenvectors GSM")



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
                                    nsigma=nsigma,
                                    number_of_points=200,
                                    distance_to_screen=100.0,
                                    scan_direction=scan_direction,
                                    sigmaxx=sigmaxx,
                                    sigmaxpxp=sigmaxpxp,
                                    useGSMapproximation=False)
    # make calculation
    out = co.calculate()

    CSD_U = out["CSD"]
    abscissas = out["abscissas"]
    eigenvalues_U = out["eigenvalues"]
    eigenvectors_U = out["eigenvectors"]
    CSD_U = out["CSD"]



    occupation_U = eigenvalues_U[0:nmodes] / eigenvalues_U.sum()
    cumulated_occupation_U = numpy.cumsum(occupation_U)
    indices = numpy.arange(abscissas.size)
    spectral_density_U = CSD_U[indices,indices]

    if do_plot:
        plot(numpy.arange(nmodes), cumulated_occupation_U, title="Cumulated occupation U")
        # plot(abscissas, spectral_density_U, title="Spectral density U")
        plot_table(abscissas, eigenvectors_U[0:10,:], title="Eigenvectors U")

    # restore spectral density from modes
    y = numpy.zeros_like(abscissas, dtype=complex)
    for i in range(200):
        y += eigenvalues_U[i] * numpy.conjugate(eigenvectors_U[i, :]) * eigenvectors_U[i, :]

    y = numpy.real(y)


    #
    # comparison
    #

    import matplotlib.pylab as plt

    ff = plot_image(CSD_GSM, 1e6*x1, 1e6*x1, title="Cross Spectral Density GSM",
                    xtitle="X1 [um]", ytitle="X2 [um]", show=0)
    plt.savefig("gsm_vs_undulator1d_figures/comparison_%s_CSD_GSM.png" % scan_direction)
    plt.show()

    ff = plot_image(numpy.abs(CSD_U), 1e6*x1, 1e6*x1, title="Cross Spectral Density U",
                    xtitle="X1 [um]", ytitle="X2 [um]", show=0)
    plt.savefig("gsm_vs_undulator1d_figures/comparison_%s_CSD_U.png" % scan_direction)
    plt.show()


    f = plot(1e6*x1, spectral_density_GSM / spectral_density_GSM.max(),
         1e6*abscissas, spectral_density_U / spectral_density_U.max(),
         xtitle="X [um]",ytitle="Spectral Density %s" % scan_direction,
         legend=["GSM", "UNDULATOR"], show=0)
    plt.savefig("gsm_vs_undulator1d_figures/comparison_%s_SD.png" % scan_direction)
    plt.show()

    f = plot(numpy.arange(nmodes), cumulated_occupation_GSM,
         numpy.arange(nmodes), cumulated_occupation_U,
         xtitle="Mode index", ytitle="Cumulated occupation",
         title="Cumulated occupation", legend=["GSM","UNDULATOR"], show=0)
    plt.savefig("gsm_vs_undulator1d_figures/comparison_%s_CumulatedO.png" % scan_direction)
    plt.show()

    for i in range(4):
        ff = plot(1e6*x1, eigenvectors_GSM[i, :] / (eigenvectors_GSM[i, :]).max(),
             1e6*abscissas, eigenvectors_U[i, :] / (eigenvectors_U[i, :]).max(),
             xtitle="X [um]", ytitle="Eigenvector",
             title="%s Mode: %d" % (scan_direction, i), legend=["GSM","UNDULATOR"], show=0)

        plt.savefig("gsm_vs_undulator1d_figures/comparison_%s_eigenvector%d.png" % (scan_direction, i))
        plt.show()

