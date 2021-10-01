
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


def run_source():
    ##########  SOURCE ##########


    #
    # create output_wavefront
    #
    #
    from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import UndulatorCoherentModeDecomposition1D
    coherent_mode_decomposition = UndulatorCoherentModeDecomposition1D(
        electron_energy=6,
        electron_current=0.2,
        undulator_period=0.018,
        undulator_nperiods=138,
        K=1.85108,
        photon_energy=7000,
        abscissas_interval=0.00025,
        number_of_points=800,
        distance_to_screen=100,
        scan_direction='V',
        sigmaxx=5.2915e-06,
        sigmaxpxp=1.88982e-06,
        useGSMapproximation=False,)
    # make calculation
    coherent_mode_decomposition_results = coherent_mode_decomposition.calculate()

    mode_index = 0
    output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(mode_index)


    if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')

    return output_wavefront

    ##########  OPTICAL SYSTEM ##########

def run_beamline(input_wavefront, f2=10.0):

    ##########  OPTICAL ELEMENT NUMBER 1 ##########



    # input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 36 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=36.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 10.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')


    #
    #---- plots -----
    #
    if plot_from_oe <= 1: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 1')


    ##########  OPTICAL ELEMENT NUMBER 2 ##########



    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape=Rectangle(-0.00025335, 0.00025335, -0.00025335, 0.00025335)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOGaussianSlit1D
    optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from_oe <= 2: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 2')


    ##########  OPTICAL ELEMENT NUMBER 3 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 29 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=29.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 1.5)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')


    #
    #---- plots -----
    #
    if plot_from_oe <= 3: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 3')


    ##########  OPTICAL ELEMENT NUMBER 4 ##########



    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.lens import WOLens1D

    optical_element = WOLens1D.create_from_keywords(
        name='',
        shape=1,
        radius=0.0005886,
        lens_aperture=0.001,
        wall_thickness=5e-05,
        material='Be',
        number_of_curved_surfaces=2,
        n_lenses=1,
        error_flag=0,
        error_file='<none>',
        error_edge_management=0,
        write_profile_flag=0,
        write_profile='profile1D.dat',
        mis_flag=0,
        xc=0,
        ang_rot=0,
        wt_offset_ffs=0,
        offset_ffs=0,
        tilt_ffs=0,
        wt_offset_bfs=0,
        offset_bfs=0,
        tilt_bfs=0)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from_oe <= 4: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 4')


    ##########  OPTICAL ELEMENT NUMBER 5 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 105 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=105.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.25)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='INTEGRAL_1D')


    #
    #---- plots -----
    #
    if plot_from_oe <= 5: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 5')


    ##########  OPTICAL ELEMENT NUMBER 6 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='',focal_length=f2)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from_oe <= 6: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 6')


    ##########  OPTICAL ELEMENT NUMBER 7 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 30 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=30.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 1.0)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='INTEGRAL_1D')


    #
    #---- plots -----
    #
    if plot_from_oe <= 7: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 7')
    return output_wavefront

if __name__ == "__main__":

    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import Tally

    if False:
        tally = Tally()

        F2 = numpy.linspace(15.0,250.0,300)
        source = run_source()

        for i,f2 in enumerate(F2):
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", i)
            output_wavefront = run_beamline(source.duplicate(), f2=f2)
            tally.append(output_wavefront)


        # tally.plot_fwhm()
        # tally.plot_intensity_at_center()
        tally.save_scan("tmp.dat")

    a = numpy.loadtxt("tmp.dat")

    plot(a[:, 0], a[:, 1], title="fwhm")
    # plot(a[:, 0], a[:, 2], title="total_intensity ")
    plot(a[:, 0], a[:, 3], title="on_axis_intensity")
    plot(a[:, 0], a[:, 4], title="peak_intensity")


