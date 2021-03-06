import numpy
from srxraylib.plot.gol import plot, plot_image, set_qt
set_qt()

def W(x1,x2):
    delta_x = x2 - x1
    mu = numpy.exp(-delta_x ** 2 / 2 / sigma_xi ** 2)
    s1 = numpy.exp(-x1 ** 2 / 2 / sigma_x ** 2)
    s2 = numpy.exp(-x2 ** 2 / 2 / sigma_x ** 2)
    return mu * numpy.sqrt(s1) * numpy.sqrt(s2)


def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF

if __name__ == "__main__":
    beta = 0.0922395 #1.151 #0.02 #
    sigma_x = 3.03783e-05

    sigma_xi = beta * sigma_x

    x1 = numpy.linspace(-0.00012, 0.00012, 400)

    # plot(x1, numpy.exp(-x1**2/(2*sigma_x**2)))

    X1 = numpy.outer(x1, numpy.ones_like(x1))
    X2 = numpy.outer(numpy.ones_like(x1), x1 )

    cross_spectral_density = W(X1,X2)
    indices = numpy.arange(x1.size)
    plot(x1, W(x1,x1),
         x1, cross_spectral_density[indices, indices],
         x1, numpy.exp(-x1**2/2/sigma_x**2),
         title="Spectral density", legend=["SD function", "SD array", "Gaussian with sigma_x"])
    plot_image(cross_spectral_density, x1, x1)

    #
    # diagonalize the CSD
    #

    w, v = numpy.linalg.eig(cross_spectral_density)
    print(w.shape, v.shape)
    idx = w.argsort()[::-1]  # large to small
    eigenvalues  = numpy.real(w[idx])
    eigenvectors = -v[:, idx].T  # minus sogn???

    print(eigenvalues[0:10])
    plot(numpy.arange(eigenvalues.size), eigenvalues)

    print("CF=", eigenvalues[0] / eigenvalues.sum(), get_coherent_fraction_exact(beta))

    from wofry.propagator.util.gaussian_schell_model import GaussianSchellModel1D, GaussianSchellModel2D


    mode = 5
    gsm = GaussianSchellModel1D(1.0, sigma_x, sigma_xi)
    y1 = eigenvectors[mode, :]
    y1 = y1 / numpy.sqrt((y1**2).sum() * (x1[1]-x1[0]))
    y2 = gsm.phi(mode, x1)
    print("modulus: ", (y2**2).sum() * (x1[1]-x1[0]))
    plot(x1, y1 ,
        x1, y2, legend=["numeric","theoretical"] )

    #
    # spectral density
    #
    sd = numpy.zeros_like(x1)
    for i in range(sd.size):
        sd[i] = cross_spectral_density[i,i]

    sdmodes = numpy.zeros_like(x1, dtype=complex)
    for i in range(sdmodes.size):
        sdmodes += eigenvalues[i] * eigenvectors[i, :]**2
    plot(x1, sd,
         x1, sdmodes,
         legend=["SD","SD from modes"],
         title="Spectral density")
