
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

def fit_profile(photon_energy):
    a = numpy.loadtxt('profile1D.dat')
    n = a.shape[0]
    w = n // 20

    x = a[(n // 2 - w):(n // 2 + w), 0] * 1e-6
    y = a[(n // 2 - w):(n // 2 + w), 1] * 1e-6

    yder = numpy.gradient(y, x)
    coeff = numpy.polyfit(x, yder, 1)

    print("\n\n\n ==========  fitted radius in the profile center : ")
    radius = 2e6 / coeff[0]
    print("fitted lens (with two curved sides) of radius = %g m " % (radius))


    import scipy.constants as codata
    import xraylib


    element = "Be"
    density = xraylib.ElementDensity(4)

    refraction_index = xraylib.Refractive_Index(element, photon_energy/1000, density)
    refraction_index_delta = 1 - refraction_index.real
    # att_coefficient = 4*numpy.pi * (xraylib.Refractive_Index(element, photon_energy/1000, density)).imag / wave_length

    f2 = radius / 2 / refraction_index_delta
    print("which corresponds for Be to a focal length of %g m " % (f2))
    return f2


#
# SOURCE========================
#




def run_source(my_mode_index=0):    
    
    
    ##########  SOURCE ##########
    
    
    #
    # create output_wavefront
    #
    #

    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-5e-05, x_max=5e-05,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(7000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=5.84299e-06, amplitude=1, mode_x=0, shift=0, beta=1.56094)

    return output_wavefront


#
# BEAMLINE========================
#




def run_beamline(output_wavefront,slit=3.62724e-05, f1=71.8241, f2=None, my_mode_index=0):
    
    
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
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=35.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
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
    
    
    ##########  OPTICAL ELEMENT NUMBER 2 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape=Rectangle(-slit/2, slit/2, -slit/2, slit/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOGaussianSlit1D
    optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)
    
    # no drift in this element 
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    
    
    ##########  OPTICAL ELEMENT NUMBER 3 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
    
    optical_element = WOScreen1D()
    
    # drift_before 31 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=31.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
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
    
    
    ##########  OPTICAL ELEMENT NUMBER 4 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    if f1 < 1000:

        optical_element = WOIdealLens1D(name='IdealLensF=71.8',focal_length=f1)

        # no drift in this element
        output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    else: # slip trace
        output_wavefront = input_wavefront

        ##########  OPTICAL ELEMENT NUMBER 5 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape=Rectangle(-0.5, 0.5, -0.5, 0.5)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)
    
    # drift_before 104 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=104.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 2.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')
    
    
    ##########  OPTICAL ELEMENT NUMBER 6 ##########
    
    
    if my_mode_index == 0: # write file only for first coherent mode
        input_wavefront = output_wavefront.duplicate()
        from orangecontrib.esrf.wofry.util.thin_object_corrector import WOThinObjectCorrector1D #TODO update

        optical_element = WOThinObjectCorrector1D(
            name='',
            file_with_thickness_mesh_flag=1,
            file_with_thickness_mesh='profile1D.dat',
            material='Be',
            focus_at=30,
            wall_thickness=0,
            apply_correction_to_wavefront=0)

        # no drift in this element
        output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    
    
    ##########  OPTICAL ELEMENT NUMBER 7 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.thin_object import WOThinObject1D #TODO update

    if f2 is None:
        f2used = fit_profile(photon_energy=input_wavefront.get_photon_energy())
    else:
        f2used = f2
    # optical_element = WOThinObject1D(name='',file_with_thickness_mesh='profile1D.dat',material='Be')
    optical_element = WOIdealLens1D(name='IdealLensF=%g' % f2used, focal_length=f2used)
    # optical_element = WOThinObject1D(name='', file_with_thickness_mesh='profile1D.dat', material='External', refraction_index_delta=6.96076e-06, att_coefficient=0)
    
    # no drift in this element 
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    
    
    ##########  OPTICAL ELEMENT NUMBER 8 ##########
    
    
    
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
    propagation_parameters.set_additional_parameters('magnification_x', 0.05)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')
    return output_wavefront


#
# MAIN FUNCTION========================
#




def main(slit=3.62724e-05, f1=71.8241, f2=None):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes
    
    tally = TallyCoherentModes()
    for my_mode_index in range(20):
        output_wavefront = run_source(my_mode_index=my_mode_index)
        output_wavefront = run_beamline(output_wavefront,slit=slit, f1=f1, f2=f2, my_mode_index=my_mode_index)
        tally.append(output_wavefront)
    
    # tally.plot_cross_spectral_density()
    # tally.plot_spectral_density()
    # tally.plot_occupation()

    return tally

#
# MAIN========================
#




# main()