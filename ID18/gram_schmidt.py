import numpy
from oasys.util.oasys_util import get_fwhm
from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from srxraylib.plot.gol import set_qt, plot_table, plot
from wofry_tools import Score

# see also https://people.inf.ethz.ch/gander/papers/qrneu.pdf


def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF

def run_wofry(xmode=0):
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=xmode, shift=0, beta=beta)

    return output_wavefront

if __name__ == "__main__":

    set_qt()

    sc = Score(scan_variable_name='mode index',additional_stored_variable_names=['mode_vector'])

    beta = 0.0922395

    for xmode in range(11):

        output_wavefront = run_wofry(xmode=xmode)

        if xmode == 0:
            WF = output_wavefront.duplicate()
        else:
            intens = WF.get_intensity()
            intens += output_wavefront.get_intensity()
            WF.set_complex_amplitude(numpy.sqrt(intens))

        intens_i = output_wavefront.get_intensity()
        # intens_i /= intens_i.sum() # normalize
        sc.append(WF, scan_variable_value=xmode, additional_stored_values=[intens_i])

    #  nagisha shirai
    # retrieve arrays
    abscissas = WF.get_abscissas()
    tmp = sc.additional_stored_values
    input_array = numpy.array(tmp)
    input_array=input_array.reshape((input_array.shape[0],input_array.shape[2]))
    normalization_factors = numpy.zeros(input_array.shape[0])

    input_array_normalized = input_array.copy()
    input_array_normalized2 = input_array.copy()
    for i in range(input_array_normalized.shape[0]):
        normalization_factors[i] = input_array_normalized[i,:].sum()
        input_array_normalized[i,:] = input_array_normalized[i,:] / input_array_normalized[i,:].sum()
        input_array_normalized2[i, :] = numpy.sqrt(input_array[i, :]) / numpy.sqrt(input_array[i, :]).sum()

    # prepare a Gaussian (data to fit)
    sigma = 0.1 * (abscissas[-1] - abscissas[0])
    # u = 0.4 * numpy.exp( - abscissas**2 / 2 / sigma**2)
    u = WF.get_intensity()
    mask = None # u

    # some tests
    # (51, 1000)
    print("Exact CF: ", get_coherent_fraction_exact(beta))
    print("From modes, CF: ", input_array.shape, normalization_factors[0] / normalization_factors.sum())


    #
    ufit1 = numpy.zeros_like(abscissas)
    projections = numpy.zeros(input_array.shape[0])
    for i in range(input_array.shape[0]):
        projections[i] = ((input_array_normalized2[i,:])**2 * (u) ).sum()
        ufit1 += normalization_factors[i] * input_array_normalized[i,:] # eigenvalue * |phi|^2
        print(i,input_array_normalized[i,:].shape , u.shape, normalization_factors[i],projections[i], normalization_factors[i]/projections[i],)
    plot(abscissas, u,
         abscissas, ufit1,
         title="test decomposition",
         legend=["data","decomposition"])

    plot(numpy.arange(input_array.shape[0]), normalization_factors,
         numpy.arange(input_array.shape[0]), (projections),
         legend=["decomposition","projections"])

    # abscissas_step = (abscissas[1] - abscissas[0])
    # for i in range(input_array_normalized.shape[0]-1):
    #     print("%d %g %g" % (
    #         i,
    #         input_array[i,:].sum() / input_array.sum(),
    #         (input_array_normalized[i,:] * input_array_normalized[i+1,:]).sum(),
    #     ))


    #
    # gram-schmidt
    #
    from orangecontrib.wofry.als.util.axo import orthonormalize_a, linear_2dgsfit1, linear_basis
    # from orangecontrib.wofry.als.util.axo_fit_profile import axo_fit_profile, calculate_orthonormal_basis

    plot_table(abscissas, numpy.sqrt(input_array_normalized), title="influence functions")


    a = []
    for i in range(input_array.shape[0]):
        # a.append({'a': input_array[i,:]/input_array[i,:].sum(), 'total_squared':0})
        # a.append({'a': input_array_normalized[i, :], 'total_squared': 0})
        a.append({'a': numpy.sqrt(input_array_normalized[i, :]), 'total_squared': 0})





    # compute the basis
    b, matrix = orthonormalize_a(a, mask=mask)

    # # normalize ????
    # for i in range(input_array.shape[1]):
    #     b[i]["a"] = b[i]["a"] / (b[i]["a"]).sum()


    # plot basis
    b_array = numpy.zeros((input_array.shape[0],input_array.shape[1]))

    for i in range(input_array.shape[0]):
        b_array[i,:] = b[i]["a"]
    plot_table(abscissas, b_array, title="basis functions")
    print(">>>> basis shape: ",b_array.shape)
    # for i in range(b_array.shape[1]-1):
    #     print(i,b_array[:,i].sum(), (b_array[:,i] * b_array[:,i+1]).sum() )


    # perform the fit
    v = linear_2dgsfit1(u, b, mask=mask)
    print("coefficients for orthonormal basis: ",v)

    vinfl = numpy.dot(matrix,v)

    print(matrix)
    print("coefficients for influence functions basis: ", vinfl.shape,vinfl)



    # evaluate the fitted data form coefficients and basis
    y = linear_basis(v, b)

    # evaluate the fitted data form coefficients and basis
    yinfl = linear_basis(vinfl, a)

    plot(abscissas, u, abscissas, y, legend=["Data", "Fit (orthonormal)"])
    plot(abscissas, u, abscissas, y, abscissas, yinfl, legend=["Data","Fit (orthonormal)","Fit (sorted influence)"])

    print(v.shape, v.max(), v.max() / v.sum())

    print(vinfl.shape, vinfl.max(), vinfl.max() / vinfl.sum())
    plot(numpy.arange(input_array.shape[0]), normalization_factors / normalization_factors.sum(),
         numpy.arange(input_array.shape[0]), v  / v.sum(),
         legend=["normalization","v"])