
#
# Import section
#
import numpy
from srxraylib.plot.gol import plot, plot_image
from score import Score



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

def run_wofry_1d(plot_from_oe=1000, mode_x=0, f1=28.2, f2=39.7):

    # plot_from_oe = 100 # set to a large number to avoid plots


    ##########  SOURCE ##########


    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012,x_max=0.00012,number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05,amplitude=1,mode_x=mode_x,shift=0,beta=0.0922395)


    if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')


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
    if plot_from_oe <= 1: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 1')


    ##########  OPTICAL ELEMENT NUMBER 2 ##########



    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    if open_slit:
        boundary_shape = Rectangle(-1.25e0, 1.25e0, -1.25e0, 1.25e0)
    else:
        boundary_shape=Rectangle(-1.25e-05, 1.25e-05, -1.25e-05, 1.25e-05)
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
    if plot_from_oe <= 3: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 3')


    ##########  OPTICAL ELEMENT NUMBER 4 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='IdealLensF=28.2',focal_length=f1)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from_oe <= 4: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 4')


    ##########  OPTICAL ELEMENT NUMBER 5 ##########



    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape=Rectangle(-0.0005, 0.0005, -0.0005, 0.0005)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # drift_before 99 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=99.0+move_f2_lense,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
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


    #
    #---- plots -----
    #
    if plot_from_oe <= 5: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 5')


    ##########  OPTICAL ELEMENT NUMBER 6 ##########
    if f2 is None: #

        input_wavefront = output_wavefront.duplicate()
        from orangecontrib.esrf.wofry.util.thin_object_corrector import WOThinObjectCorrector1D  # TODO update

        optical_element = WOThinObjectCorrector1D(
            name='',
            file_with_thickness_mesh_flag=1,
            file_with_thickness_mesh='profile1D.dat',
            material='Be',
            focus_at=36-move_f2_lense,
            wall_thickness=0,
            apply_correction_to_wavefront=1)

        # no drift in this element
        output_wavefront = optical_element.applyOpticalElement(input_wavefront)
        # F2.append("?")
        f22 = 0

    elif f2 == 0:

        input_wavefront = output_wavefront.duplicate()
        from orangecontrib.esrf.wofry.util.thin_object_corrector import WOThinObjectCorrector1D  # TODO update

        optical_element = WOThinObjectCorrector1D(
            name='',
            file_with_thickness_mesh_flag=1,
            file_with_thickness_mesh='profile1D.dat',
            material='Be',
            focus_at=36-move_f2_lense,
            wall_thickness=0,
            apply_correction_to_wavefront=0)

        # no drift in this element
        output_wavefront = optical_element.applyOpticalElement(input_wavefront)

        f22 = fit_profile(output_wavefront)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> f22: ", f22)
        # F2.append(f22)


        #  now apply the ideal lens with the calculated f2
        from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

        optical_element = WOIdealLens1D(name='IdealLensF=39.7',focal_length=f22)

        # no drift in this element
        output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    else:

        input_wavefront = output_wavefront.duplicate()
        from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

        optical_element = WOIdealLens1D(name='IdealLensF=39.7',focal_length=f2)

        # no drift in this element
        output_wavefront = optical_element.applyOpticalElement(input_wavefront)
        # F2.append(f2)
        f22 = f2
    #
    #---- plots -----
    #
    if plot_from_oe <= 6: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 6')


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
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=36.000000 - move_f2_lense,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    if move_f2_lense > 0:
        propagation_parameters.set_additional_parameters('magnification_x', 0.05)
    else:
        propagation_parameters.set_additional_parameters('magnification_x', 0.25)

    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')


    #
    #---- plots -----
    #
    if plot_from_oe <= 7: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 7')

    return output_wavefront, f22


def run_multimode(up_to_mode=0, f1=28.2, f2=None):
    for i in range(up_to_mode+1):

        if i == 0:
            wf, tmp_f2 = run_wofry_1d(plot_from_oe=1000, mode_x=i, f1=f1, f2=f2)
            WF = wf.duplicate()
            F2 = tmp_f2
        else:
            if f2 is None:
                wf, tmp_f2 = run_wofry_1d(plot_from_oe=1000, mode_x=i, f1=f1, f2=F2)
            else:
                if f2 == 0:
                    wf, tmp_f2 = run_wofry_1d(plot_from_oe=1000, mode_x=i, f1=f1, f2=F2)
                else:
                    wf, tmp_f2 = run_wofry_1d(plot_from_oe=1000, mode_x=i, f1=f1, f2=f2)

            intens = WF.get_intensity()
            intens += wf.get_intensity()
            WF.set_complex_amplitude(numpy.sqrt(intens))
    return WF, F2


def get_f2(f1=28.2, verbose=False):
    position_source = 0
    position_slit = 35.0
    position_lens1 = 65.0
    position_lens2 = 164.0
    position_sample = 200.0

    p1 = position_lens1 - position_source
    q1 = 1 / (1 / f1 - 1 / p1)
    D = position_lens2 - position_lens1
    p2 = D - q1

    q2 = position_sample - position_lens2
    f2 = 1.0 / (1 / p2 + 1 / q2)

    if verbose:
        print("\n\nD: %g, q1+p2: %g" % (D, q1 + p2))
        print("p1: %g" % p1)
        print("q1: %g" % q1)
        print("p2: %g" % p2)
        print("q2: %g" % q2)
        print("Total length: %g" % (p1 + q1 + p2 + q2))

    return f2


def get_refraction_index(photon_energy=10000.0):
    import scipy.constants as codata
    import xraylib
    wave_length = codata.h * codata.c / codata.e / photon_energy

    element = "Be"
    density = xraylib.ElementDensity(4)

    refraction_index = xraylib.Refractive_Index(element, photon_energy/1000, density)
    refraction_index_delta = 1 - refraction_index.real
    att_coefficient = 4*numpy.pi * (xraylib.Refractive_Index(element, photon_energy/1000, density)).imag / wave_length

    return refraction_index_delta, att_coefficient


def fit_profile(wf):
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
    photon_energy = wf.get_photon_energy()

    element = "Be"
    density = xraylib.ElementDensity(4)

    refraction_index = xraylib.Refractive_Index(element, photon_energy/1000, density)
    refraction_index_delta = 1 - refraction_index.real
    # att_coefficient = 4*numpy.pi * (xraylib.Refractive_Index(element, photon_energy/1000, density)).imag / wave_length

    f2 = radius / 2 / refraction_index_delta
    print("which corresponds for Be to a focal length of %g m " % (f2))
    return f2

if __name__ == "__main__":

    do_plot = True
    save_file = True
    up_to_mode = 0
    open_slit = 1

    move_f2_lense = 24.0
    f2 = 0    # this means guess it and use the ideal lens.
    # f2 = 39.7 # this is fixed   29.38 #
    # f2 = None # this means use the corrector and do not compute ideal lens


    #
    #
    #

    sc = Score(scan_variable_name='f1 [m]',additional_stored_variable_names=['f2 [m]'])


    F1 = numpy.linspace(18.2, 38.2, 5)
    # F1 = numpy.linspace(28.2,28.3,5)
    F1 = numpy.linspace(10, 40, 100)
    # F2 = []




    if do_plot:
            WF1, F2_1 = run_multimode(up_to_mode=up_to_mode, f1=F1[0], f2=f2)
            WF2, F2_2 = run_multimode(up_to_mode=up_to_mode, f1=F1[F1.size // 2], f2=f2)
            WF3, F2_3 = run_multimode(up_to_mode=up_to_mode, f1=F1[-1], f2=f2)

            if f2 is None:
                ff2_1 = "f2=?"
                ff2_2 = "f2=?"
                ff2_3 = "f2=?"
            else:
                ff2_1 = "f2=%g" % F2_1
                ff2_2 = "f2=%g" % F2_2
                ff2_3 = "f2=%g" % F2_3
            plot(WF1.get_abscissas(), WF1.get_intensity(),
                 WF2.get_abscissas(), WF2.get_intensity(),
                 WF3.get_abscissas(), WF3.get_intensity(),
                 legend=["f1=%g %s"%(F1[0],F2_1), "f1=%g %s"%(F1[F1.size//2],F2_2), "f1=%g %s"%(F1[-1],F2_3)],
                 title="")


    for i, f1 in enumerate(F1):
        WF, F2 = run_multimode(up_to_mode=up_to_mode, f1=f1, f2=f2)

        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> F2, f2", F2, f2)
        sc.append(WF, scan_variable_value=f1, additional_stored_values=[F2])


    if open_slit:
        filename = "tmp_openslit_uptomode%d.dat"% up_to_mode
    else:
        filename = "tmp_uptomode%d.dat" % up_to_mode

    sc.save(filename)


    sc.plot()

    source_fwhm = 15.14e-6
    plot(numpy.array(sc.scan_variable_value),
         numpy.array(sc.fwhm) / source_fwhm,
         title="Magnification")

    # plot(numpy.array(sc.scan_variable_value),
    #      numpy.array(F2),
    #      title="F1 F2 map")

    # f = open("tmp_f1f2.dat",'w')
    # for i in range(len(F2)):
    #     f.write("%g %g\n" % (sc.scan_variable_value[i], F2[i]))
    # f.close()
