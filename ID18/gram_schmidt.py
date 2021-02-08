import numpy
from oasys.util.oasys_util import get_fwhm
from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from srxraylib.plot.gol import set_qt, plot_table, plot, plot_image
from wofry_tools import Score

# see also https://people.inf.ethz.ch/gander/papers/qrneu.pdf


def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF

def run_wofry_source(xmode=0):
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=xmode, shift=0, beta=beta)

    return output_wavefront

def run_wofry_to_screen(xmode=0, slit_size = 25e-6):
    #
    # Import section
    #
    import numpy

    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

    from wofryimpl.propagator.propagators1D.fresnel import Fresnel1D
    from wofryimpl.propagator.propagators1D.fresnel_convolution import FresnelConvolution1D
    from wofryimpl.propagator.propagators1D.fraunhofer import Fraunhofer1D
    from wofryimpl.propagator.propagators1D.integral import Integral1D
    from wofryimpl.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
    from wofryimpl.propagator.propagators1D.fresnel_zoom_scaling_theorem import FresnelZoomScaling1D

    from srxraylib.plot.gol import plot, plot_image
    plot_from_oe = 100  # set to a large number to avoid plots

    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=xmode, shift=0, beta=0.0922395)

    if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

    ##########  OPTICAL SYSTEM ##########

    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 35 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=35.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 8.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    #
    # ---- plots -----
    #
    if plot_from_oe <= 1: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 1')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle

    boundary_shape = Rectangle(-slit_size/2, slit_size/2, -slit_size/2, slit_size/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    #
    # ---- plots -----
    #
    if plot_from_oe <= 2: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 2')

    return output_wavefront
if __name__ == "__main__":

    set_qt()

    sc = Score(scan_variable_name='mode index',additional_stored_variable_names=['mode_vector'])

    beta = 0.0922395

    for xmode in range(51):

        output_wavefront = run_wofry_source(xmode=xmode)
        # output_wavefront = run_wofry_to_screen(xmode=xmode, slit_size=1.0)

        if xmode == 0:
            WF = output_wavefront.duplicate()
        else:
            intens = WF.get_intensity()
            intens += output_wavefront.get_intensity()
            WF.set_complex_amplitude(numpy.sqrt(intens))

        complex_amplitude_i = output_wavefront.get_complex_amplitude()
        sc.append(WF, scan_variable_value=xmode, additional_stored_values=[complex_amplitude_i])

    #  nagisha shirai
    # retrieve arrays
    abscissas = WF.get_abscissas()
    tmp = sc.additional_stored_values
    input_array = numpy.array(tmp)
    input_array=input_array.reshape((input_array.shape[0],input_array.shape[2]))

    eigenvalues_old = numpy.zeros(input_array.shape[0])
    eigenvectors_old  = numpy.zeros_like(input_array, dtype=float)
    for i in range(eigenvalues_old.size):
        eigenvalues_old[i] = ((input_array[i, :]) ** 2).sum()
        eigenvectors_old[i, :] = numpy.real(input_array[i, :] / numpy.sqrt(eigenvalues_old[i]))

    # check orthonormality, CF
    for i in range(eigenvalues_old.size - 1):
        print(i,
              (eigenvectors_old[i, :] * eigenvectors_old[i, :]).sum(),
              (eigenvectors_old[i, :] * eigenvectors_old[i + 1, :]).sum(),
              (eigenvectors_old[i, :] * eigenvectors_old[i, :] * eigenvectors_old[i + 1, :] * eigenvectors_old[i + 1, :]).sum(),
              )
    print("Exact CF: ", get_coherent_fraction_exact(beta))
    print("From modes, CF: ", eigenvalues_old[0] / eigenvalues_old.sum())


    # prepare a Gaussian (data to fit)
    sigma = 0.1 * (abscissas[-1] - abscissas[0])
    # u = 0.4 * numpy.exp( - abscissas**2 / 2 / sigma**2)
    u = WF.get_intensity()
    mask = None # u

    # # some tests
    # ufit1 = numpy.zeros_like(abscissas)
    # projections = numpy.zeros(input_array.shape[0])
    # I = numpy.zeros(input_array.shape[0])
    #
    # for i in range(eigenvalues_old.size):
    #     # projections[i] = ((eigenvectors[i,:])**2 * (u) ).sum()
    #     projections[i] = ((eigenvectors_old[i, :]) * (u)).sum()
    #     I[i] = (eigenvectors_old[i, :]).sum()
    #     ufit1 += eigenvalues_old[i] * (numpy.real(eigenvectors_old[i, :])) ** 2
    #     # print(i, eigenvalues[i], projections[i])
    #
    # plot(abscissas, u,
    #      abscissas, ufit1,
    #      title="test decomposition",
    #      legend=["data","decomposition"])
    #
    # try1 = projections * I
    # for i in range(eigenvalues_old.size):
    #     print(">>", i, eigenvalues_old[i], projections[i], try1[i])
    #
    # plot(numpy.arange(eigenvalues_old.size), eigenvalues_old,
    #      # numpy.arange(eigenvalues.size), projections ** 2 / projections[0] ** 2,
    #      numpy.arange(eigenvalues_old.size), try1,
    #      legend=["eigenvalues","try1"])
    #
    # print(">>>>", u.sum(), eigenvalues_old.sum(), try1.sum())

    cross_spectral_density = numpy.zeros((abscissas.size, abscissas.size))


    for i in range(eigenvalues_old.size):
        cross_spectral_density += eigenvalues_old[i] * numpy.outer(eigenvectors_old[i, :], eigenvectors_old[i, :])

    plot_image(cross_spectral_density)
    print("matrix shape: ", cross_spectral_density.shape)
    w, v = numpy.linalg.eig(cross_spectral_density)
    print(w.shape, v.shape)
    idx = w.argsort()[::-1]  # large to small
    eigenvalues  = numpy.real(w[idx])
    eigenvectors = numpy.real(v[:, idx].T)
    for i in range(10):
        print(eigenvalues_old[i], eigenvalues[i])
    plot(abscissas, eigenvectors_old[0, :],
         abscissas, eigenvectors[0, :], )

    print("Exact CF: ", get_coherent_fraction_exact(beta))
    print("From modes, CF_old: ", eigenvalues_old[0] / eigenvalues_old.sum())
    print("From modes, CF    : ", eigenvalues[0] / eigenvalues.sum())

    # abscissas_step = (abscissas[1] - abscissas[0])
    # for i in range(input_array_normalized.shape[0]-1):
    #     print("%d %g %g" % (
    #         i,
    #         input_array[i,:].sum() / input_array.sum(),
    #         (input_array_normalized[i,:] * input_array_normalized[i+1,:]).sum(),
    #     ))


    # #
    # # gram-schmidt
    # #
    # from orangecontrib.wofry.als.util.axo import orthonormalize_a, linear_2dgsfit1, linear_basis
    # # from orangecontrib.wofry.als.util.axo_fit_profile import axo_fit_profile, calculate_orthonormal_basis
    #
    # plot_table(abscissas, numpy.sqrt(input_array_normalized), title="influence functions")
    #
    #
    # a = []
    # for i in range(input_array.shape[0]):
    #     # a.append({'a': input_array[i,:]/input_array[i,:].sum(), 'total_squared':0})
    #     # a.append({'a': input_array_normalized[i, :], 'total_squared': 0})
    #     a.append({'a': numpy.sqrt(input_array_normalized[i, :]), 'total_squared': 0})
    #
    #
    #
    #
    #
    # # compute the basis
    # b, matrix = orthonormalize_a(a, mask=mask)
    #
    # # # normalize ????
    # # for i in range(input_array.shape[1]):
    # #     b[i]["a"] = b[i]["a"] / (b[i]["a"]).sum()
    #
    #
    # # plot basis
    # b_array = numpy.zeros((input_array.shape[0],input_array.shape[1]))
    #
    # for i in range(input_array.shape[0]):
    #     b_array[i,:] = b[i]["a"]
    # plot_table(abscissas, b_array, title="basis functions")
    # print(">>>> basis shape: ",b_array.shape)
    # # for i in range(b_array.shape[1]-1):
    # #     print(i,b_array[:,i].sum(), (b_array[:,i] * b_array[:,i+1]).sum() )
    #
    #
    # # perform the fit
    # v = linear_2dgsfit1(u, b, mask=mask)
    # print("coefficients for orthonormal basis: ",v)
    #
    # vinfl = numpy.dot(matrix,v)
    #
    # print(matrix)
    # print("coefficients for influence functions basis: ", vinfl.shape,vinfl)
    #
    #
    #
    # # evaluate the fitted data form coefficients and basis
    # y = linear_basis(v, b)
    #
    # # evaluate the fitted data form coefficients and basis
    # yinfl = linear_basis(vinfl, a)
    #
    # plot(abscissas, u, abscissas, y, legend=["Data", "Fit (orthonormal)"])
    # plot(abscissas, u, abscissas, y, abscissas, yinfl, legend=["Data","Fit (orthonormal)","Fit (sorted influence)"])
    #
    # print(v.shape, v.max(), v.max() / v.sum())
    #
    # print(vinfl.shape, vinfl.max(), vinfl.max() / vinfl.sum())
    # plot(numpy.arange(input_array.shape[0]), normalization_factors / normalization_factors.sum(),
    #      numpy.arange(input_array.shape[0]), v  / v.sum(),
    #      legend=["normalization","v"])