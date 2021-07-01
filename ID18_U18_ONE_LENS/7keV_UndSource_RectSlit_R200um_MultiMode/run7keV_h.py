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


def run_source(my_mode_index=0,number_of_points=800, zero_emittance=0):
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

        if zero_emittance:
            sigmaxx=1e-07
            sigmaxpxp=1e-08
        else:
            sigmaxx=2.97321e-05
            sigmaxpxp=4.37237e-06

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
            sigmaxx=sigmaxx,
            sigmaxpxp=sigmaxpxp,
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


def run_beamline(output_wavefront,aperture=40e-6, distance=30.0):
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
    propagation_parameters.set_additional_parameters('magnification_x', 8.0)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='INTEGRAL_1D')

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
    propagation_parameters.set_additional_parameters('magnification_x', 2.5)
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

    ##########  OPTICAL ELEMENT NUMBER 5 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-0.5, 0.5, -0.5, 0.5)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # drift_before 18.4 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=distance,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.1)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='INTEGRAL_1D')

    return output_wavefront


#
# MAIN FUNCTION========================
#


def main(aperture=40e-6, distance=30.0, number_of_points=800,nmodes=10):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()


    for my_mode_index in range(nmodes):
        print(">>>>>>>>>>>>>>>>>>>>> mode %d of %d" % (my_mode_index, nmodes))
        output_wavefront = run_source(my_mode_index=my_mode_index,number_of_points=number_of_points)
        output_wavefront = run_beamline(output_wavefront, aperture=aperture, distance=distance)
        tally.append(output_wavefront)

    # tally.plot_cross_spectral_density(show=1, filename="")
    # tally.plot_spectral_density(show=1, filename="")
    # tally.plot_occupation(show=1, filename="")

    return tally


#
# MAIN========================
#


# main()
if __name__ == "__main__":
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes, Tally
    from oasys.util.oasys_util import get_fwhm
    from srxraylib.plot.gol import plot
    #
    #
    #
    # size_at_aperture = 565e-6
    APERTURE = [40.3e-6, 85.1e-6, 145e-6, 1000e-6, -40.3e-6, -85.1e-6, -145e-6, -1000e-6] # [ 5000e-6] # [-40.3e-6, -85.1e-6, -145e-6, -1000e-6] #
    DISTANCE = numpy.linspace(10, 50, 50) # numpy.array([18.4]) #   # # 31.19 28.4
    number_of_points = 800 # 800


    for aperture in APERTURE:

        # src1, wf = main(aperture=aperture, distance=18.4168, number_of_points=number_of_points)

        filename = "aperture_h_%g.dat" % (1e6 * aperture) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        f = open(filename, 'w')

        f.write("# S 1 scored data\n")
        f.write("# N 5\n")
        f.write("# L  distance  fwhm  total_intensity  on_axis_intensity  peak_intensity")

        if aperture < 0:
            aperture *= -1
            nmodes = 1
        else:
            nmodes = 10

        for i,distance in enumerate(DISTANCE):
            tally = main(aperture=aperture, distance=distance, nmodes=nmodes)

            spectral_density = tally.get_spectral_density() # numpy.zeros_like(abscissas)
            abscissas = tally.get_abscissas()
            fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)

            I = spectral_density
            x = abscissas

            fwhm, quote, coordinates = get_fwhm(I, x)
            intensity_at_center = I[I.size // 2]
            intensity_total = I.sum() * (x[1] - x[0])
            intensity_peak = I.max()


            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # plot(1e6 * abscissas, spectral_density,
            #      legend=["From Cross Spectral Density"],
            #      xtitle="x [um]", ytitle="Spectral Density", title="D=%g m,FWHM = %g um, a=%g um" % (distance, fwhm, aperture*1e6))

            f.write("\n %g  %g  %g  %g  %g  " % (distance,  fwhm,  intensity_total,  intensity_at_center,  intensity_peak))

        f.close()
        print("File %s written to disk" % filename)
        # tally.save("aperture_h_%g.dat" % (aperture))


    # main()
