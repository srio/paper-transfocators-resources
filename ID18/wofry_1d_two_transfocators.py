import numpy

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
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(7000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.04613e-05, amplitude=1, mode_x=0, shift=0, beta=0.118608)
    # previous command is useless but...
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.04613e-05, amplitude=1, mode_x=my_mode_index, shift=0,
                                               beta=0.118608)
    return output_wavefront


#
# BEAMLINE========================
#


def run_beamline(output_wavefront,
                f2 = 143.648000,
                f1 = 71.824100):
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
    propagation_parameters.set_additional_parameters('magnification_x', 4.0)
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
    boundary_shape = Rectangle(-1.81362e-05, 1.81362e-05, -1.81362e-05, 1.81362e-05)
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
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=31.000000, q=0.000000,
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
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D


    optical_element = WOIdealLens1D(name='IdealLensF=71.8', focal_length=f1)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    ##########  OPTICAL ELEMENT NUMBER 5 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-0.0005, 0.0005, -0.0005, 0.0005)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # drift_before 104 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=104.000000, q=0.000000,
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

    ##########  OPTICAL ELEMENT NUMBER 6 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D


    optical_element = WOIdealLens1D(name='IdealLensF=143.648', focal_length=f2)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    # ##########  OPTICAL ELEMENT NUMBER 7 ##########
    #
    # input_wavefront = output_wavefront.duplicate()
    # from orangecontrib.esrf.wofry.util.thin_object_corrector import WOThinObjectCorrector1D  # TODO update
    #
    # optical_element = WOThinObjectCorrector1D(
    #     name='',
    #     file_with_thickness_mesh_flag=1,
    #     file_with_thickness_mesh='profile1D.dat',
    #     material='Be',
    #     focus_at=30,
    #     wall_thickness=0,
    #     apply_correction_to_wavefront=1)
    #
    # # no drift in this element
    # output_wavefront = optical_element.applyOpticalElement(input_wavefront)

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

def run_all():
    #
    # MAIN========================
    #


    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(50):
        output_wavefront = run_source(my_mode_index=my_mode_index)
        output_wavefront = run_beamline(output_wavefront)
        tally.append(output_wavefront)

    cf, eigenvalues, eigenvectors, cross_spectral_density = tally.calculate_coherent_fraction(do_plot=True)
    print('Coherent fraction from new (rediagonalized) modes: %f ' % cf)


if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    directions = ['h','v']
    focii = ['large', 'small', ]
    energies = [7, 10, 15, 20, 35]
    distances = [170, 192]

    import datastorage
    data = datastorage.read("summary_of_GSM_results_for_Manuel.npz")


    for distance in distances:
        for energy in energies:
            for focus in focii:
                label = "e%02dkeV_f2_at_%dm_%s" % (energy, distance, focus)
                for direction in directions:
                    # print("%s  %s" % (label, direction))
                    tmp = data[label][direction]
                    for key in tmp.keys():
                        # print(key, tmp[key])
                        slit      = tmp["p035m_hard_aperture"]
                        f1        = tmp["p066m_f1_used"]
                        f2        = tmp["pf2pos_f2_used"]
                        size200   = tmp["p200m_fwhm_beam"]
                        sizewaist = tmp["waist_fwhm_size"]
                        poswaist  = tmp["waist_pos"]
                    print(">>>>%s %s :  slit: %g m, f1: %g, f2: %g, size: %g um (%g at %g m)" %
                          (label, direction, slit, f1, f2, 1e6*size200, 1e6*sizewaist, poswaist))



    # attached the npz file, best is to read it with "datastorage": data = datastorage.read(fname)
    #
    #
    # and here the δ,β
    #
    # In [1]: from sr import materials
    #
    # In [2]: for energy in (7,10,15,20,35):
    #    ...: print(materials.get_delta_beta("Be",energy=energy))
    # ...:
    # (6.959562521724472e-06, 3.920823829597519e-09)
    # (3.4079264780162433e-06, 1.1184654031792269e-09)
    # (1.514042063277543e-06, 3.6034468029877136e-10)
    # (8.515233095307551e-07, 2.0080453378588767e-10)
    # (2.780114856104632e-07, 8.816014655748295e-11) """



    # print(data["e07keV_f2_at_170m_large"]['v'])

    # F1 = []
    # F2 = []
    # for label in key1:
    #     # label = "e07keV_f2_at_170m_large"
    #     for direction in directions:
    #         tmp = data[label][direction]
    #         for key in tmp.keys():
    #             # print(key, tmp[key])
    #             slit      = tmp["p035m_hard_aperture"]
    #             f1        = tmp["p066m_f1_used"]
    #             f2        = tmp["pf2pos_f2_used"]
    #             size200   = tmp["p200m_fwhm_beam"]
    #             sizewaist = tmp["waist_fwhm_size"]
    #             poswaist  = tmp["waist_pos"]
    #         F1.append(f1)
    #         F2.append(f2)
    #         print(">>>>%s %s :  slit: %g m, f1: %g, f2: %g, size: %g um (%g at %g m)" %
    #               (label, direction, slit, f1, f2, 1e6*size200, 1e6*sizewaist, poswaist))

    # print(ENERGY)
    # print(F1)
    # print(len(ENERGY),len(F1))
    # plot(numpy.array(F1),numpy.array(F2)) #,xrange=[0,100], yrange=[0,100])
    # print(key1[0], data[key1[0]]['h'])

    # run_all()





