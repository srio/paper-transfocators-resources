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

        # output_wavefront = run_wofry_source(xmode=xmode)
        output_wavefront = run_wofry_to_screen(xmode=xmode, slit_size=25e-6)


        complex_amplitude_i = output_wavefront.get_complex_amplitude()

        if xmode == 0:
            WF = output_wavefront.duplicate()
            intens = output_wavefront.get_intensity()
            WF.set_complex_amplitude(numpy.sqrt(intens))
        else:
            intens = WF.get_intensity()
            intens += output_wavefront.get_intensity()
            WF.set_complex_amplitude(numpy.sqrt(intens))

        sc.append(WF, scan_variable_value=xmode, additional_stored_values=[complex_amplitude_i])


    # retrieve arrays
    abscissas = WF.get_abscissas()
    tmp = sc.additional_stored_values
    input_array = numpy.zeros((len(tmp), abscissas.size), dtype=complex)
    for i in range(len(tmp)):
        input_array[i,:] = tmp[i][0]
    nmodes = input_array.shape[0]

    #
    # calculate and diagonalize the CSD
    #
    cross_spectral_density = numpy.zeros((abscissas.size, abscissas.size), dtype=complex)


    for i in range(nmodes):
        cross_spectral_density += numpy.outer(numpy.conjugate(input_array[i, :]), input_array[i, :])

    plot_image(numpy.abs(cross_spectral_density))
    print("matrix cross_spectral_density: ", cross_spectral_density.shape)
    w, v = numpy.linalg.eig(cross_spectral_density)
    print(w.shape, v.shape)
    idx = w.argsort()[::-1]  # large to small
    eigenvalues  = numpy.real(w[idx])
    eigenvectors = v[:, idx].T


    print("Exact CF at the source: ", get_coherent_fraction_exact(beta))
    print("From modes, CF    : ", eigenvalues[0] / eigenvalues.sum())

    #
    # plot intensity
    #
    abscissas_step = (abscissas[1] - abscissas[0])
    y = numpy.zeros_like(abscissas)

    for i in range(nmodes):
        y += eigenvalues[i] * numpy.real(numpy.conjugate(eigenvectors[i,:]) * eigenvectors[i,:])


    plot(abscissas, WF.get_intensity(),
         abscissas, y, legend=["Data", "From modes"])

    plot(numpy.arange(nmodes), eigenvalues[0:nmodes]/(eigenvalues[0:nmodes].sum()),title="mode occupation")
