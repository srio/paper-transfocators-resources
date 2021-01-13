#
#
#

import numpy
from srxraylib.plot.gol import plot


from wofry.propagator.util.gaussian_schell_model import GaussianSchellModel1D, GaussianSchellModel2D


def get_eigenvalue(beta,index):

    A = 1.0
    sigma_i = 1.0
    sigma_mu = beta * sigma_i

    aa = 1.0 / (4 * sigma_i**2)
    bb = 1.0 / (2 * sigma_mu**2)
    cc = numpy.sqrt( aa**2 + 2*aa*bb)

    return A * numpy.sqrt(numpy.pi/(aa+bb+cc)) * (bb/(aa+bb+cc))**index
#
def get_eigenvalues_sum(beta,index_max):
    out = 0.0
    for i in range(index_max+1):
        out += get_eigenvalue(beta,i)
    return out
#
# def get_coherent_fraction(beta,index_max):
#     return get_eigenvalue(beta,0) / get_eigenvalues_sum(beta,index_max)
#
# def get_coherent_fraction_exact(beta):
#     q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
#     q = 1.0 / q
#     CF = 1 - q
#     return CF
#
def get_eigenvalues_array(beta,index_max):
    modes = numpy.arange(index_max)
    lambdai = numpy.zeros(index_max)
    for i in range(lambdai.size):
        lambdai[i] = get_eigenvalue(beta,i) # A * numpy.sqrt(numpy.pi/(aa+bb+cc)) * (bb/(aa+bb+cc))**i
    occupation = lambdai/get_eigenvalues_sum(beta,index_max-1)
    cumulated_occupation = occupation.cumsum()

    return modes,lambdai,occupation,cumulated_occupation

def sum_all_eigenvalues_over_lambda0(beta):
    q = 1.0 / (1 + beta**2 / 2 + beta * numpy.sqrt(1 + (beta / 2)**2))
    return 1.0 / (1.0 - q)

def sum_n_eigenvalues_over_lambda0(beta, n):
    q = 1.0 / (1 + beta**2 / 2 + beta * numpy.sqrt(1 + (beta / 2)**2))
    return (1.0 - q**n) / (1.0 - q)

def cumulated_occupation_n_eigenvalues(beta, n):
    q = 1.0 / (1 + beta ** 2 / 2 + beta * numpy.sqrt(1 + (beta / 2) ** 2))
    return (1.0 - q**n)

def number_of_modes(beta, cumulated_occupation=0.9):
    q = 1.0 / (1 + beta ** 2 / 2 + beta * numpy.sqrt(1 + (beta / 2) ** 2))
    # print(">>>>", cumulated_occupation, 1.0 - cumulated_occupation)
    return int( numpy.log(1.0 - cumulated_occupation) / numpy.log(q))

def occupation(beta, n):
    q = 1.0 / (1 + beta ** 2 / 2 + beta * numpy.sqrt(1 + (beta / 2) ** 2))
    return q**n * ( 1 - q)


    # return sum_n_eigenvalues_over_lambda0(beta, n) / sum_all_eigenvalues_over_lambda0(beta)

# def sortedModeIndices(eigenvalues_x, eigenvalues_y, index_energy, n_points=50):
#     # if self._sorted_mode_indices is None:
#     #     eigenvalues_x = np.array([self._mode_x.beta(i) for i in range(n_points)])
#     #     eigenvalues_y = np.array([self._mode_y.beta(i) for i in range(n_points)])
#     #
#     #     f = np.outer(eigenvalues_x, eigenvalues_y)
#     #     self._sorted_mode_indices = np.array(np.unravel_index(f.flatten().argsort()[::-1], (n_points, n_points)))
#
#     f = numpy.outer(eigenvalues_x, eigenvalues_y)
#     _sorted_mode_indices = numpy.array(numpy.unravel_index(f.flatten().argsort()[::-1], (n_points, n_points)))
#
#     n, m = _sorted_mode_indices[:, index_energy]
#     return n, m

if __name__ == "__main__":

# ---- Vartanyants, photon energy=10000.000000 eV
# H sigmaI= 30.38, sigmaMu=  2.83, beta=  0.09, CF_ratio=  0.09, CF_GSM=  0.09:
# V sigmaI=  5.00, sigmaMu=  3.76, beta=  0.75, CF_ratio=  0.67, CF_GSM=  0.52:
# ---- srio, photon energy=10000.000000 eV
# H sigmaI= 30.38, sigmaMu=  2.80, beta=  0.09, CF_ratio=  0.09, CF_GSM=  0.09:
# V sigmaI=  5.00, sigmaMu=  5.76, beta=  1.15, CF_ratio=  0.67, CF_GSM=  0.67:

    beta_h = 0.09
    beta_v = 1.15
    sigma_s_h = 30.38e-6
    sigma_s_v = 5.00e-6

    # method 0
    modes_h, lambdai, occupation_h, cumulated_occupation_h = get_eigenvalues_array(beta_h,99)
    modes_v, lambdai, occupation_v, cumulated_occupation_v = get_eigenvalues_array(beta_v,99)


    # f = plot(
    #     modes_h, occupation_h,
    #     modes_v, occupation_v,xtitle="mode index",ytitle="occupation",
    #     legend=["H beta=%4.1f" % beta_h, "V beta=%4.1f" % beta_v],
    #     xrange=[0,20],show=True,figsize=(10,7.5),marker=['o','o'])

    plot(
        modes_h, cumulated_occupation_h,
        modes_v, cumulated_occupation_v, xtitle="mode index", ytitle="cumulated occupation",
        legend=["H beta=%4.1f" % beta_h, "V beta=%4.1f" % beta_v],
        xrange=[0, 100], show=True)

    # method 1
    gsm_h = GaussianSchellModel1D(1.0, sigma_s_h, sigma_s_h * beta_h)
    gsm_v = GaussianSchellModel1D(1.0, sigma_s_v, sigma_s_h * beta_v)

    modes_h = numpy.arange(100)
    for i in range(10):
        print("i, lambdai (mg), lambdai (here): ", i, gsm_h.beta(i), gsm_h.beta(0) * (1 + beta_h**2 / 2 + beta_h * numpy.sqrt(1 + (beta_h/2)**2))**(-i))

    print("\n\n\n\n\n")
    lambda0 = gsm_h.beta(0)
    sum_eigenvalues = 0.0
    for i in range(10):
        sum_eigenvalues += gsm_h.beta(i)
        print("i, sum eigenvalues eq, array, previous: ", i,
              1e6 * lambda0 * sum_n_eigenvalues_over_lambda0(beta_h, i+1),
              1e6 * sum_eigenvalues,
              1e6 * cumulated_occupation_h[i] * (lambda0 * sum_all_eigenvalues_over_lambda0(beta_h)))

    print("\n\n\n\n\n")

    for i in range(90):

        print("i, cumulated old, cumulated_new: ", i,
              cumulated_occupation_h[i],
              cumulated_occupation_n_eigenvalues(beta_h, i+1),
              number_of_modes(beta_h, cumulated_occupation_h[i]))


    for i in range(10):
        print("occupation: ", i, occupation_h[i], occupation(beta_h, i))

    #
    # sort indices
    #
    limit_cumulated_occupation = 0.99
    n_h = number_of_modes(beta_h, limit_cumulated_occupation)
    n_v = number_of_modes(beta_v, limit_cumulated_occupation)
    print("Needed number of harmonics H, V: ",n_h, n_v)

    gsm_2d = GaussianSchellModel2D(1.0, sigma_s_h, sigma_s_h * beta_h, sigma_s_v, sigma_s_v * beta_v)
    sum_eigenvalues = 0.0
    up_to_i = n_h * n_v


    sorted_array = []
    for i in range(up_to_i):
        ih, iv = gsm_2d.sortedModeIndices(i, n_points=up_to_i)
        # ch = cumulated_occupation_n_eigenvalues(beta_h,ih+1)
        # cv = cumulated_occupation_n_eigenvalues(beta_v,iv+1)
        sum_eigenvalues += gsm_2d.beta(ih, iv)
        sorted_array.append([ih, iv])


        # eigenvalues_x = numpy.array([gsm_2d._mode_x.beta(i) for i in range(up_to_i)])
        # eigenvalues_y = numpy.array([gsm_2d._mode_y.beta(i) for i in range(up_to_i)])
        # f = numpy.outer(eigenvalues_x, eigenvalues_y)
        # _sorted_mode_indices = numpy.array(numpy.unravel_index(f.flatten().argsort()[::-1], (up_to_i, up_to_i)))
        # ih, iv = _sorted_mode_indices[:, i]

        print(i, ih, iv, 1e12 * sum_eigenvalues)


    sorted_array2 = numpy.array(sorted_array)
    print("Maximum harmonics indices used: ", sorted_array2[:,0].max(), sorted_array2[:,1].max() )




