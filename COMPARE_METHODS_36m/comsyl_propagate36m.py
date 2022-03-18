
#
# Import section
#
import numpy

from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D

from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
from wofryimpl.propagator.propagators2D.fresnel import Fresnel2D
from wofryimpl.propagator.propagators2D.fresnel_convolution import FresnelConvolution2D
from wofryimpl.propagator.propagators2D.fraunhofer import Fraunhofer2D
from wofryimpl.propagator.propagators2D.integral import Integral2D
from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D

from srxraylib.plot.gol import plot, plot_image
plot_from_oe = 100 # set to a large number to avoid plots


def get_result_for_single_mode(af_oasys, mode_index=0):
    ##########  SOURCE ##########


    #
    # create output_wavefront
    #
    #
    # from comsyl.autocorrelation.CompactAFReader import CompactAFReader
    # filename = "/scisoft/users/srio/mycomsyl/calculations/comsyl_EBS_ID18_CPM18_7keV_s3.0.npz"
    # af_oasys = CompactAFReader.initialize_from_file(filename)

    # mode_index = 0
    output_wavefront = af_oasys.get_wavefront(mode_index,normalize_with_eigenvalue=1)


    if plot_from_oe <= 0: plot_image(output_wavefront.get_intensity(),output_wavefront.get_coordinate_x(),output_wavefront.get_coordinate_y(),aspect='auto',title='SOURCE')


    ##########  OPTICAL SYSTEM ##########





    ##########  OPTICAL ELEMENT NUMBER 1 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen

    optical_element = WOScreen()

    # drift_before 36 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=0.5*2.5+36.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
    propagation_parameters.set_additional_parameters('magnification_x', 2 * 4.0)
    propagation_parameters.set_additional_parameters('magnification_y', 3 * 10.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoomXY2D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_XY_2D')

    #
    # ---- plots -----
    #

    if plot_from_oe <= 1: plot_image(output_wavefront.get_intensity(), output_wavefront.get_coordinate_x(),
                                     output_wavefront.get_coordinate_y(), aspect='auto', title='OPTICAL ELEMENT NR 1')

    return output_wavefront

def add_modes():
    from comsyl.autocorrelation.CompactAFReader import CompactAFReader
    filename = "/scisoft/users/srio/mycomsyl/calculations/comsyl_EBS_ID18_CPM18_7keV_s3.0.npz"
    af_oasys = CompactAFReader.initialize_from_file(filename)

    for i in range(250):
        output_wavefront = get_result_for_single_mode(af_oasys, mode_index=i)

        if i == 0:
            INTENSITY = output_wavefront.get_intensity()
        else:
            INTENSITY += output_wavefront.get_intensity()

    return INTENSITY, output_wavefront.get_coordinate_x(), output_wavefront.get_coordinate_y()

#
#================ MAIN ===================
#

if __name__ == "__main__":


    INTENSITY, x, y = add_modes()

    plot_image(INTENSITY, 1e6*x, 1e6*y, aspect='auto',title='mode')