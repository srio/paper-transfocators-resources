
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

from wofry_tools import Score
#
#
#
from syned.beamline.shape import *

def run_wofry_1d(plot_from=0, mode_x=0):
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012,x_max=0.00012,number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05,amplitude=1,mode_x=mode_x,shift=0,beta=0.0922395)


    if plot_from <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')

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
    propagation_parameters.set_additional_parameters('magnification_x', 8.0)
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
    if plot_from <= 1:  plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 1')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########



    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *

    # boundary_shape=Rectangle(-1.25e-05, 1.25e-05, -1.25e-05, 1.25e-05)
    boundary_shape = Rectangle(-1e-6*slit_size_in_um/2, 1e-6*slit_size_in_um/2, -1e-6*slit_size_in_um/2, 1e-6*slit_size_in_um/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D, WOGaussianSlit1D
    if use_gaussian_slits:
        optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)
    else:
        optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from <= 2: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 2')

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
    if plot_from <= 3: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 3')

    ##########  OPTICAL ELEMENT NUMBER 4 ##########



    input_wavefront = output_wavefront.duplicate()


    if use_real_lens:
        from orangecontrib.esrf.wofry.util.lens import WOLens1D
        optical_element = WOLens1D.create_from_keywords(
            name='',
            shape=1,
            radius=0.000192435,
            lens_aperture=0.001,
            wall_thickness=5e-05,
            material='Be',
            refraction_index_delta=5.3e-07, # used if material='External'
            att_coefficient=0.00357382, # used if material='External'
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
    else:
        from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D
        optical_element = WOIdealLens1D(name='IdealLensF=28.2',focal_length=28.200000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from <= 4: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 4')

    ##########  OPTICAL ELEMENT NUMBER 5 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_after 99 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=0.000000,    q=q5,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #

    propagation_parameters.set_additional_parameters('magnification_x', magnification_x)
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
    if plot_from <= 5: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 5')

    return output_wavefront

def run_multimode(up_to_mode=0):
    for i in range(up_to_mode+1):
        wf = run_wofry_1d(plot_from=1000, mode_x=i)
        if i == 0:
            WF = wf.duplicate()
        else:
            intens = WF.get_intensity()
            intens += wf.get_intensity()
            WF.set_complex_amplitude(numpy.sqrt(intens))
    return WF

if __name__ == "__main__":


    use_gaussian_slits = False
    use_real_lens = False
    npoints = 100
    do_plot = True
    save_file = False
    up_to_mode = 0

    NSIGMAS = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 4.0, 6.0]


    for nsigmas in NSIGMAS:
        sc = Score(scan_variable_name='q5 [m]')

        if nsigmas <= 2:
            magnification_x = 2
        if nsigmas <= .5:
            magnification_x = 1
        else:
            magnification_x = 4

        slit_size_in_um = 125.0 / 2.35 * nsigmas
        Q5 = numpy.concatenate((
            numpy.linspace(10, 70, npoints // 3),
            numpy.linspace(71, 300, npoints // 3),  # 49.8098 at source, 470 at slit
            numpy.linspace(301, 500, npoints // 3) ))

        npoints = Q5.size


        if do_plot:
            q5 = Q5[0] + 0.1
            WF1 = run_multimode(up_to_mode=up_to_mode)
            q5 = 49.8098
            WF2 = run_multimode(up_to_mode=up_to_mode)
            q5 = 470
            WF3 = run_multimode(up_to_mode=up_to_mode)
            q5 = Q5[-1]
            WF4 = run_multimode(up_to_mode=up_to_mode)

            plot(WF1.get_abscissas(), WF1.get_intensity(),
                 WF2.get_abscissas(), WF2.get_intensity(),
                 WF3.get_abscissas(), WF3.get_intensity(),
                 WF4.get_abscissas(), WF4.get_intensity(),
                 legend=["%g"%Q5[0], "%g"%50, "%g"%470, "%g"%Q5[-1]],
                 title="nsigma = %g" % nsigmas)

        for i in range(npoints):
            # run_wofry_1d(plot_from=50)
            q5 = Q5[i]

            WF = run_multimode(up_to_mode=up_to_mode)
            sc.append(WF, scan_variable_value=q5)

            if numpy.mod(i,10) == 0:
                print(">>>>>> iteration index %d of %d" % (i,npoints+1))
                # if 0:
                #     plot(WF.get_abscissas(), WF.get_intensity(),
                #          title=">>>>>> iteration index %d of %d" % (i,npoints))

        if use_real_lens:
            use_real_lens_extension = "R"
        else:
            use_real_lens_extension = ""

        if use_gaussian_slits:
            use_gaussian_slits_extension = "G"
        else:
            use_gaussian_slits_extension = ""

        if up_to_mode > 0:
            up_to_mode_extension = "M"
        else:
            up_to_mode_extension = ""

        if save_file:
            sc.save(filename="tmp%s%s%s%2.1f.dat" % (up_to_mode_extension,use_gaussian_slits_extension,use_real_lens_extension,nsigmas))

        if do_plot: sc.plot()