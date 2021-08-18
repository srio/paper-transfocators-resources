
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
plot_from_oe = 100 # set to a large number to avoid plots


def run_beamline(ang_rot=0.1, distance=29.3087, ):
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.001, x_max=0.001,
                                                                          number_of_points=2000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_plane_wave_from_complex_amplitude(complex_amplitude=complex(1, 0), inclination=0)

    if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

    ##########  OPTICAL SYSTEM ##########

    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    aperture = 0.001
    boundary_shape = Rectangle(-aperture/2, aperture/2, -aperture/2, aperture/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    #
    # ---- plots -----
    #
    if plot_from_oe <= 1: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 1')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.lens import WOLens1D

    optical_element = WOLens1D.create_from_keywords(
        name='',
        shape=1,
        radius=0.0005,
        lens_aperture=0.005,
        wall_thickness=1e-05,
        material='Be',
        number_of_curved_surfaces=2,
        n_lenses=1,
        error_flag=0,
        error_file='<none>',
        error_edge_management=0,
        write_profile_flag=0,
        write_profile='profile1D.dat',
        mis_flag=1,
        xc=0,
        ang_rot=ang_rot,
        wt_offset_ffs=0,
        offset_ffs=0,
        tilt_ffs=0,
        wt_offset_bfs=0,
        offset_bfs=0,
        tilt_bfs=0)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    #
    # ---- plots -----
    #
    if plot_from_oe <= 2: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 2')

    ##########  OPTICAL ELEMENT NUMBER 3 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 73.2716 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=distance,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.2)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='INTEGRAL_1D')

    #
    # ---- plots -----
    #
    if plot_from_oe <= 3: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 3')


    return output_wavefront

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image, set_qt
    from orangecontrib.esrf.wofry.util.tally import Tally
    set_qt()

    f0 = 73.2716
    Factor = numpy.linspace(0.5, 1.3, 55)  # 1.0 - 0.3
    tally = Tally(do_store_wavefronts=True)
    ang_rot = 0.3

    for factor in Factor:
        output_wavefront = run_beamline(ang_rot=ang_rot, distance=f0*factor)
        tally.append(output_wavefront,scan_variable_value=factor)

    min_index = tally.get_fwhm().argmin()
    min_size = 1e6 * tally.get_fwhm()[min_index] #.min()
    min_position = f0 * Factor[min_index]

    tally.plot_fwhm(factor_fwhm=1e6,
                    title="tilt=%3.2f rad; Minimum size: %g um, at d=%g m, d/f0=%g " % \
                                                (ang_rot, min_size, min_position, Factor[min_index]),
                    xtitle="distance/f0",
                    )
    tally.plot_wavefronts_intensity(factor_abscissas=1e6,
                                    xtitle="distance/f0",
                                    ytitle="wavefront abscissas [um]",
                                    title="Intensity for tilt = %3.2f rad; f0=%g m" % (ang_rot, f0))






