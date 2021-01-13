#
#
#

import numpy
from srxraylib.plot.gol import plot


def get_eigenvalue(beta,index):

    A = 1.0
    sigma_i = 1.0
    sigma_mu = beta * sigma_i


    aa = 1.0 / (4 * sigma_i**2)
    bb = 1.0 / (2 * sigma_mu**2)
    cc = numpy.sqrt( aa**2 + 2*aa*bb)

    return A * numpy.sqrt(numpy.pi/(aa+bb+cc)) * (bb/(aa+bb+cc))**index

def get_eigenvalues_sum(beta,index_max):
    out = 0.0
    for i in range(index_max+1):
        out += get_eigenvalue(beta,i)
    return out

def get_coherent_fraction(beta,index_max):
    return get_eigenvalue(beta,0) / get_eigenvalues_sum(beta,index_max)

def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF

def get_eigenvalues_array(beta,index_max):
    modes = numpy.arange(index_max)
    lambdai = numpy.zeros(index_max)
    for i in range(lambdai.size):
        lambdai[i] = get_eigenvalue(beta,i) # A * numpy.sqrt(numpy.pi/(aa+bb+cc)) * (bb/(aa+bb+cc))**i
    occupation = lambdai/get_eigenvalues_sum(beta,index_max-1)
    cumulated_occupation = occupation.cumsum()

    return modes,lambdai,occupation,cumulated_occupation


if __name__ == "__main__":

    modes, lambdai, occupation_0p1, cumulated_occupation_0p1 = get_eigenvalues_array(0.1,99)
    modes, lambdai, occupation_1,   cumulated_occupation_1   = get_eigenvalues_array(1.0,99)
    modes, lambdai, occupation_10,  cumulated_occupation_10  = get_eigenvalues_array(10.0,99)

    f = plot(modes, occupation_0p1,
        modes, occupation_1,
        modes, occupation_10,xtitle="mode index",ytitle="occupation",
        legend=["beta=0.1","beta=1.0","beta=10.0"],
        xrange=[0,25],show=True,figsize=(10,7.5),marker=['o','o','o'])


    plot(modes, cumulated_occupation_0p1,
        modes, cumulated_occupation_1,
        modes, cumulated_occupation_10,xtitle="mode index",ytitle="cumulated occupation",
        legend=["beta=0.1, CF=%f"%cumulated_occupation_0p1[0],
                "beta=1.0, CF=%f"%cumulated_occupation_1[0],
                "beta=10.0, CF=%f"%cumulated_occupation_10[0]],
        xrange=[0,25],show=True)


    # import matplotlib.pylab as plt
    # plt.bar(modes, cumulated_occupation_1, linewidth=0.1)
    # plt.bar(modes, cumulated_occupation_0p1, linewidth=0.8)
    # plt.show()
    #
    # coherent fraction
    #

    betas = numpy.linspace(0.,10.0, 100)
    CF1 = numpy.zeros_like(betas)
    CF = numpy.zeros_like(betas)

    for i,beta in enumerate(betas):
        CF[i] = get_coherent_fraction_exact(beta)
        CF1[i] = get_coherent_fraction(beta,1000)
        print(i,CF[i],CF1[i])

    f = plot(betas,CF,xtitle=r'$\mathrm{\beta}$',ytitle="CF",show=True,figsize=(10,7.5))

    # outfile = "cf_gaussianshell.eps"
    # plt.savefig(outfile)
    # print("File written to disk: %s"%outfile)

    # plt.show()

