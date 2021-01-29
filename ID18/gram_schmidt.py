if __name__ == "__main__":
    import numpy
    from oasys.util.oasys_util import get_fwhm
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
    from srxraylib.plot.gol import set_qt, plot_table, plot
    from wofry_tools import Score
    set_qt()

    sc = Score(scan_variable_name='mode index',additional_stored_variable_names=['mode_vector'])


    for xmode in range(11):
        output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                              number_of_points=1000)
        output_wavefront.set_photon_energy(10000)
        output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=xmode, shift=0, beta=0.0922395)

        if xmode == 0:
            WF = output_wavefront.duplicate()
        else:
            intens = WF.get_intensity()
            intens += output_wavefront.get_intensity()
            WF.set_complex_amplitude(numpy.sqrt(intens))

        intens_i = output_wavefront.get_intensity()
        intens_i /= intens_i.sum() # normalize
        sc.append(WF, scan_variable_value=xmode, additional_stored_values=[intens_i])


    #
    # gram-schmidt
    #
    from orangecontrib.wofry.als.util.axo import orthonormalize_a, linear_2dgsfit1, linear_basis
    # from orangecontrib.wofry.als.util.axo_fit_profile import axo_fit_profile, calculate_orthonormal_basis

    tmp = sc.additional_stored_values
    input_array = numpy.array(tmp)
    input_array=input_array.reshape((input_array.shape[0],input_array.shape[2]))
    input_array = input_array.T
    abscissas=WF.get_abscissas()
    plot_table(abscissas, input_array.T)


    a = []
    for i in range(input_array.shape[1]):
        a.append({'a': input_array[:, i], 'total_squared':0})


    # prepare a Gaussian (data to fit)
    sigma = 0.1 * (abscissas[-1] - abscissas[0])
    u = 0.4 * numpy.exp( - abscissas**2 / 2 / sigma**2)

    mask = None # u


    # compute the basis
    b, matrix = orthonormalize_a(a, mask=mask)


    # plot basis
    b_array = numpy.zeros((input_array.shape[0],input_array.shape[1]))

    for i in range(input_array.shape[1]):
        b_array[:,i] = b[i]["a"]
    plot_table(abscissas, b_array.T, title="basis functions")


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

    plot(abscissas, u, abscissas, y,legend=["Data", "Fit (orthonormal)"])
    # plot(abscissas,u,abscissas,y,abscissas,yinfl,legend=["Data","Fit (orthonormal)","Fit (sorted influence)"])