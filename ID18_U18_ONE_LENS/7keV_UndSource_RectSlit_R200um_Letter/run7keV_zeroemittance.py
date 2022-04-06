import numpy

def run_wofry(source=None, aperture=1.0, distance=18.4168, number_of_points=2000):
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
    plot_from_oe = 50  # set to a large number to avoid plots

    ##########  SOURCE ##########

    if source is None:
        print(">>>>>>>>>>>>>>>>>> CALCULATING...")
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
            number_of_points=number_of_points,
            distance_to_screen=100,
            scan_direction='V',
            sigmaxx=1e-07,
            sigmaxpxp=1e-08,
            useGSMapproximation=False, )
        # make calculation
        coherent_mode_decomposition_results = coherent_mode_decomposition.calculate()

        mode_index = 0
        output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(mode_index)
    else:
        print(">>>>>>>>>>>>>>>>> RE-USING")
        output_wavefront = source

    ret_source = output_wavefront.duplicate()
    if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

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
    propagation_parameters.set_additional_parameters('magnification_x', 15.0)
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
    if plot_from_oe <= 1: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 1')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-aperture/2, aperture/2, -aperture/2, aperture/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

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

    #
    # ---- plots -----
    #
    if plot_from_oe <= 3: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 3')

    ##########  OPTICAL ELEMENT NUMBER 4 ##########

    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.lens import WOLens1D

    optical_element = WOLens1D.create_from_keywords(
        name='',
        shape=1,
        radius=0.0002,
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
    # ---- plots -----
    #
    if plot_from_oe <= 4: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 4')

    ##########  OPTICAL ELEMENT NUMBER 5 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 18.4168 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=distance, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.05)
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
    if plot_from_oe <= 5: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 5')

    return ret_source, output_wavefront

if __name__ == "__main__":
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes, Tally



    size_at_aperture = 565e-6
    APERTURE_FACTOR = [0.0001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 4, 6]
    DISTANCE = numpy.linspace(10, 50, 1000)
    number_of_points = 3000

    src1, wf = run_wofry(source=None, aperture=APERTURE_FACTOR[0], distance=18.4168, number_of_points=number_of_points)

    for aperture_factor in APERTURE_FACTOR:
        aperture = size_at_aperture * aperture_factor

        tally = Tally(scan_variable_name='distance')

        for i,distance in enumerate(DISTANCE):
            print(src1)
            src2, wf = run_wofry(source=src1, aperture=aperture, distance=distance)
            tally.append(wf, scan_variable_value=distance)

        # tally.plot()
        tally.save_scan("zeroemittance/aperture_factor_%g.dat" % (aperture_factor))