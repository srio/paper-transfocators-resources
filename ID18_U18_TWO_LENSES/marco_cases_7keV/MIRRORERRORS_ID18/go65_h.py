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
            scan_direction='H',
            sigmaxx=2.97321e-05,
            sigmaxpxp=4.37237e-06,
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


def run_beamline(output_wavefront):
    ##########  OPTICAL SYSTEM ##########

    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    input_wavefront = output_wavefront.duplicate()

    from orangecontrib.esrf.wofry.util.mirror import WOMirror1D

    optical_element = WOMirror1D.create_from_keywords(
        name='',
        shape=0,
        flip=0,
        p_focus=1,
        q_focus=1,
        grazing_angle_in=0.003,
        p_distance=29.5,
        q_distance=0.25,
        zoom_factor=3,
        error_flag=1,
        error_file='dabam_profile_140325822323120.dat',
        error_file_oversampling_factor=100,
        mirror_length=0,
        mirror_points=0,
        write_profile=0)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.toolbox import WOToolbox1D  # TODO update

    optical_element = WOToolbox1D(name='', crop_factor=1.2, abscissas_factor=1, shift_center=0, change_photon_energy=0,
                                  new_photon_energy=0)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    ##########  OPTICAL ELEMENT NUMBER 3 ##########

    input_wavefront = output_wavefront.duplicate()

    from orangecontrib.esrf.wofry.util.mirror import WOMirror1D

    optical_element = WOMirror1D.create_from_keywords(
        name='',
        shape=0,
        flip=1,
        p_focus=1,
        q_focus=1,
        grazing_angle_in=0.003,
        p_distance=0.5,
        q_distance=6,
        zoom_factor=1.6,
        error_flag=1,
        error_file='dabam_profile_140325822345040.dat',
        error_file_oversampling_factor=200,
        mirror_length=0,
        mirror_points=0,
        write_profile=0)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    ##########  OPTICAL ELEMENT NUMBER 4 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-2e-05, 2e-05, -2e-05, 2e-05)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    ##########  OPTICAL ELEMENT NUMBER 5 ##########

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
    propagation_parameters.set_additional_parameters('magnification_x', 2.5)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')
    return output_wavefront


#
# MAIN FUNCTION========================
#


def main():
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(30):
        output_wavefront = run_source(my_mode_index=my_mode_index)
        output_wavefront = run_beamline(output_wavefront)
        tally.append(output_wavefront)

    tally.plot_cross_spectral_density(show=1, filename="")
    tally.plot_spectral_density(show=1, filename="")
    tally.plot_occupation(show=1, filename="")

    tally.save_spectral_density(filename="tmp65_h_spectral_density.dat")
    tally.save_occupation(filename="tmp65_h_occupation.dat")


#
# MAIN========================
#


main()