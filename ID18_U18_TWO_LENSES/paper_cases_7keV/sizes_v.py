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


#
# SOURCE========================
#


def run_source(my_mode_index=0):
    global coherent_mode_decomposition
    try:
        if my_mode_index == 0: raise Exception()
        tmp = coherent_mode_decomposition
    except:

        ##########  SOURCE ##########

        #
        # create output_wavefront
        #
        #
        from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import \
            UndulatorCoherentModeDecomposition1D
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
            useGSMapproximation=False, )
        # make calculation
        coherent_mode_decomposition_results = coherent_mode_decomposition.calculate()

        mode_index = 0
        output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(mode_index)
    output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(my_mode_index)
    return output_wavefront


#
# BEAMLINE========================
#


def run_beamline(output_wavefront,aperture=40e-6, radius1=10.0, radius2=10.0):
    ##########  OPTICAL SYSTEM ##########

    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 36 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=36.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 10.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-aperture/2, aperture/2, -aperture/2, aperture/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

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
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=29.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 1.5)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    ##########  OPTICAL ELEMENT NUMBER 4 ##########

    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.lens import WOLens1D

    optical_element = WOLens1D.create_from_keywords(
        name='',
        shape=1,
        radius=radius1, #0.0006921,
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
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=105.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.5)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='INTEGRAL_1D')

    ##########  OPTICAL ELEMENT NUMBER 6 ##########

    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.lens import WOLens1D

    optical_element = WOLens1D.create_from_keywords(
        name='',
        shape=1,
        radius=radius2, #0.0005529,
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
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=30.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.25)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')
    return output_wavefront


if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image, set_qt
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes
    import xraylib
    import os

    set_qt()

    APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]
    aperture = APERTURE[0]
    direction = 'v'

    fileroot = "case7keV_%s" % direction
    subdirectory = "./sizes_slit%g_%s" % (1e6 * aperture, direction)

    os.mkdir(subdirectory)

    a = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (1e6 * aperture, direction))
    F1 = a[:, 0].copy()
    F2 = a[:, 1].copy()


    R1 = []
    R2 = []
    for i in range(F1.size):
        xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7,1.85)).real
        R1.append(F1[i] * (2 * xrl_delta))
        R2.append(F2[i] * (2 * xrl_delta))

    for i in range(len(R1)):
        tally = TallyCoherentModes()
        for my_mode_index in range(10):
            output_wavefront = run_source(my_mode_index=my_mode_index)
            output_wavefront = run_beamline(output_wavefront, aperture=aperture, radius1=R1[i], radius2=R2[i])
            tally.append(output_wavefront)

        tally.plot_spectral_density(show=0,filename="%s/%s_spectral_density_%03d.png" % (subdirectory, fileroot, i),
                                    title="R1: %d um, R2: %d um, F1: %4.1f, F2: %4.1f " % (1e6*R1[i], 1e6*R2[i], F1[i], F2[i]))
        tally.save_spectral_density(filename="%s/%s_spectral_density_%03d.dat" % (subdirectory, fileroot,i))
        tally.save_occupation(filename="%s/%s_occupation_%03d.dat" % (subdirectory, fileroot, i))