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
from oasys.util.oasys_util import get_fwhm

from syned.beamline.shape import *


def get_wavefront_intensity_fwhm(wf):
    fwhm, quote, coordinates = get_fwhm(wf.get_intensity(), wf.get_abscissas())
    return fwhm

def get_wavefront_intensity_I0(wf):
    I = wf.get_intensity()
    return I[I.size // 2]

def get_wavefront_intensity_ITOTAL(wf):
    I = wf.get_intensity()
    x = wf.get_abscissas()
    return I.sum() * (x[1] - x[0])

def run_wofry_v(plot_from=0,run_up_to_element=100,q7=38.000000,mode_x=0):
    tkt = {}
    tkt["DISTANCE"] = []
    tkt["FWHM"] = []
    tkt["I0"] = []
    tkt["ITOTAL"] = []


    ##########  SOURCE ##########


    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-5e-05,x_max=5e-05,number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=5.00125e-06,amplitude=1,mode_x=mode_x,shift=0,beta=1.15083)


    if plot_from <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(0.0)
    tkt["output_wavefront"] = output_wavefront

    if run_up_to_element == 0: return tkt


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
    propagation_parameters.set_additional_parameters('magnification_x', 16.0)
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
    if plot_from <= 1: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 1')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(35)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 1: return tkt

    ##########  OPTICAL ELEMENT NUMBER 2 ##########



    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *
    boundary_shape=Rectangle(-6.35e-05, 6.35e-05, -6.35e-05, 6.35e-05)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from <= 2: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 2')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(35.01)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 2: return tkt

    ##########  OPTICAL ELEMENT NUMBER 3 ##########



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
    propagation_parameters.set_additional_parameters('magnification_x', 0.5)
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
    if plot_from <= 3: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 3')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(55)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 3: return tkt

    ##########  OPTICAL ELEMENT NUMBER 4 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='IdealLensF=48.9',focal_length=48.900000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from <= 4: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 4')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(55.01)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 4: return tkt

    ##########  OPTICAL ELEMENT NUMBER 5 ##########



    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *
    boundary_shape=Rectangle(-0.0005, 0.0005, -0.0005, 0.0005)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # drift_before 99 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=99.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 1.0)
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
    if plot_from <= 5: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 5')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(164)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 5: return tkt

    ##########  OPTICAL ELEMENT NUMBER 6 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='IdealLensF=38.3',focal_length=38.300000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from <= 6: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 6')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(164.01)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 6: return tkt

    ##########  OPTICAL ELEMENT NUMBER 7 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 36 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=0.0,    q=q7,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.5)
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
    if plot_from <= 7: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 7')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(164.0 + q7)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 7: return tkt

    return tkt

def run_wofry_h(plot_from=0,run_up_to_element=100,q7=38.000000,mode_x=0):
    tkt = {}
    tkt["DISTANCE"] = []
    tkt["FWHM"] = []
    tkt["I0"] = []
    tkt["ITOTAL"] = []
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=mode_x, shift=0, beta=0.0922395)

    if plot_from <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(0.0)
    tkt["output_wavefront"] = output_wavefront

    if run_up_to_element == 0: return tkt
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
    if plot_from <= 1: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 1')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(35)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 1: return tkt

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *
    boundary_shape = Rectangle(-1.25e-05, 1.25e-05, -1.25e-05, 1.25e-05)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    #
    # ---- plots -----
    #
    if plot_from <= 2: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 2')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(35.01)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 2: return tkt


    ##########  OPTICAL ELEMENT NUMBER 3 ##########

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

    #
    # ---- plots -----
    #
    if plot_from <= 3: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 3')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(55)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 3: return tkt

    ##########  OPTICAL ELEMENT NUMBER 4 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='IdealLensF=28.2', focal_length=28.200000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    #
    # ---- plots -----
    #
    if plot_from <= 4: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 4')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(55.01)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 4: return tkt

    ##########  OPTICAL ELEMENT NUMBER 5 ##########

    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *
    boundary_shape = Rectangle(-0.0005, 0.0005, -0.0005, 0.0005)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # drift_before 99 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=99.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 2.0)
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
    if plot_from <= 5: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 5')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(164)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 5: return tkt

    ##########  OPTICAL ELEMENT NUMBER 6 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='IdealLensF=39.7', focal_length=39.700000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    #
    # ---- plots -----
    #
    if plot_from <= 6: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 6')
    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(164.01)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 6: return tkt

    ##########  OPTICAL ELEMENT NUMBER 7 ##########

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
                                       coordinates=ElementCoordinates(p=0.0, q=q7,
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

    #
    # ---- plots -----
    #
    if plot_from <= 7: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 7')

    tkt["FWHM"].append(get_wavefront_intensity_fwhm(output_wavefront))
    tkt["I0"].append(get_wavefront_intensity_I0(output_wavefront))
    tkt["ITOTAL"].append(get_wavefront_intensity_ITOTAL(output_wavefront))
    tkt["DISTANCE"].append(164.0 + q7)
    tkt["output_wavefront"] = output_wavefront
    if run_up_to_element == 7: return tkt

    return tkt

def run_wofry_multimode(run_up_to_mode=10,nelements=7,q7=38.000000,horizontal_or_vertical=0):

    # first run to define outputs


    FWHM     = []
    I0       = []
    ITOTAL   = []
    DISTANCE = []

    for ielement in range(nelements+1):
        for imode in range(run_up_to_mode+1):
            # print(">>>>>>> run element %d  mode %d" % (ielement, imode))
            if horizontal_or_vertical == 0:
                tkti = run_wofry_h(plot_from=1000, run_up_to_element=ielement, q7=q7, mode_x=imode)
            else:
                tkti = run_wofry_v(plot_from=1000, run_up_to_element=ielement, q7=q7, mode_x=imode)

            if imode == 0:
                abscissas = tkti["output_wavefront"].get_abscissas()
                intensity = tkti["output_wavefront"].get_intensity()
            else:
                intensity += tkti["output_wavefront"].get_intensity()

        DISTANCE.append( tkti["DISTANCE"][-1] )
        I0.append( intensity[intensity.size // 2] )
        ITOTAL.append(intensity.sum() * (abscissas[1] - abscissas[0]))
        fwhm, _, _ = get_fwhm(intensity, abscissas)
        FWHM.append( fwhm )
        # plot(abscissas, intensity, title="element %d  mode %d" % (ielement, imode))

    return DISTANCE, FWHM, I0, ITOTAL

if __name__ == "__main__":

    q7 = 36.000000
    run_up_to_element = 7

    horizontal_or_vertical = 1

    run_mode = 1 # 0 = single mode, 1=all modes

    if horizontal_or_vertical == 0:
        pscan = numpy.linspace(175, 250.0, 10)
        title = "HORIZONTAL"
        run_up_to_mode = 50
    else:
        pscan = numpy.linspace(175, 250.0, 30)
        title = "VERTICAL"
        run_up_to_mode = 5

    #
    # single mode
    #
    if run_mode == 0:
        if horizontal_or_vertical == 0:
            tkt = run_wofry_h(plot_from=7, q7=q7, run_up_to_element=run_up_to_element, mode_x=0)
        else:
            tkt = run_wofry_v(plot_from=7, q7=q7, run_up_to_element=run_up_to_element, mode_x=0)

        print(tkt["DISTANCE"])
        print(tkt["FWHM"])
        print(tkt["I0"])
        print(tkt["ITOTAL"])

        # remove last elemet
        tkt["FWHM"].pop()
        tkt["I0"].pop()
        tkt["ITOTAL"].pop()
        tkt["DISTANCE"].pop()

        FWHM = tkt["FWHM"]
        I0 = tkt["I0"]
        ITOTAL = tkt["ITOTAL"]
        DISTANCE = tkt["DISTANCE"]

        for i in range(pscan.size):
            if horizontal_or_vertical == 0:
                tkt = run_wofry_h(plot_from=18,q7=pscan[i]-164,run_up_to_element=run_up_to_element, mode_x=0)
            else:
                tkt = run_wofry_v(plot_from=18, q7=pscan[i] - 164, run_up_to_element=run_up_to_element, mode_x=0)

            DISTANCE.append(tkt["DISTANCE"][-1])
            FWHM.append(tkt["FWHM"][-1])
            I0.append(tkt["I0"][-1])
            ITOTAL.append(tkt["ITOTAL"][-1])
    elif run_mode == 1: # many modes

        DISTANCE, FWHM, I0, ITOTAL = run_wofry_multimode(run_up_to_mode=run_up_to_mode,
                                                 nelements=run_up_to_element,
                                                 q7=q7,
                                                 horizontal_or_vertical=horizontal_or_vertical)
        print(">>>> DISTANCE: ", DISTANCE)
        print(">>>> FWHM: ", FWHM)

        # remove last point (image) to prepare scan
        DISTANCE.pop()
        FWHM.pop()
        I0.pop()
        ITOTAL.pop()

        for i in range(pscan.size):
            q7i = pscan[i] - 164
            print(">>>>>>>  screen distance: %g (%d of %d)" % (q7i,i,pscan.size))
            DISTANCEi, FWHMi, I0i, ITOTALi = run_wofry_multimode(run_up_to_mode=run_up_to_mode,
                                                        nelements=run_up_to_element,
                                                        q7=q7i,
                                                        horizontal_or_vertical=horizontal_or_vertical)

            DISTANCE.append(DISTANCEi[-1])
            FWHM.append(FWHMi[-1])
            I0.append(I0i[-1])
            ITOTAL.append(ITOTALi[-1])

    f = open("tmp.dat",'w')
    for i in range(len(DISTANCE)):
        f.write("%g %g %g %g\n" % (DISTANCE[i], FWHM[i], I0[i], ITOTAL[i]))
    f.close()
    print("File written to disk: tmp.dat")


    plot(DISTANCE, I0, title=title, ytitle="Peak intensity [a.u.]", xtitle="Distance from source [m]", figsize=(15,4), show=0)
    plot(DISTANCE, ITOTAL, title=title, ytitle="Beam intensity [a.u.]", xtitle="Distance from source [m]", figsize=(15, 4), show=0)
    plot(DISTANCE, 1e6*numpy.array(FWHM), title=title, ytitle="FWHM [um]", xtitle="Distance from source [m]", figsize=(15,4), show=1)