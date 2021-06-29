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

def run_source_und_h(my_mode_index=0, energy_in_keV=7):
    global coherent_mode_decomposition_h
    try:
        tmp = coherent_mode_decomposition_h
    except:


        ##########  SOURCE ##########

        if energy_in_keV == 7:
            K = 1.85108
        elif energy_in_keV == 15:
            K = 0.729628
        elif energy_in_keV == 30:
            K = 1.341095
        else:
            raise Exception("Please provide K value")
        #
        # create output_wavefront
        #
        #
        from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import UndulatorCoherentModeDecomposition1D
        coherent_mode_decomposition_h = UndulatorCoherentModeDecomposition1D(
            electron_energy=6,
            electron_current=0.2,
            undulator_period=0.018,
            undulator_nperiods=138,
            K=K,
            photon_energy=1e3 * energy_in_keV,
            abscissas_interval=0.00025,
            number_of_points=800,
            distance_to_screen=100,
            scan_direction='H',
            sigmaxx=2.97321e-05,
            sigmaxpxp=4.37237e-06,
            useGSMapproximation=False,)
        # make calculation
        coherent_mode_decomposition_results = coherent_mode_decomposition_h.calculate()

        mode_index = 0
        output_wavefront = coherent_mode_decomposition_h.get_eigenvector_wavefront(mode_index)
    output_wavefront = coherent_mode_decomposition_h.get_eigenvector_wavefront(my_mode_index)
    return output_wavefront

def run_source_und_v(my_mode_index=0, energy_in_keV=7):
    global coherent_mode_decomposition_v
    try:
        tmp = coherent_mode_decomposition_v
    except:

        if energy_in_keV == 7:
            K = 1.85108
        elif energy_in_keV == 15:
            K = 0.729628
        elif energy_in_keV == 30:
            K = 1.341095
        else:
            raise Exception("Please provide K value")
        ##########  SOURCE ##########

        #
        # create output_wavefront
        #
        #
        from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import \
            UndulatorCoherentModeDecomposition1D
        coherent_mode_decomposition_v = UndulatorCoherentModeDecomposition1D(
            electron_energy=6,
            electron_current=0.2,
            undulator_period=0.018,
            undulator_nperiods=138,
            K=K,
            photon_energy=1e3 * energy_in_keV,
            abscissas_interval=0.00025,
            number_of_points=800,
            distance_to_screen=100,
            scan_direction='V',
            sigmaxx=5.2915e-06,
            sigmaxpxp=1.88982e-06,
            useGSMapproximation=False, )
        # make calculation
        coherent_mode_decomposition_results = coherent_mode_decomposition_v.calculate()

        mode_index = 0
        output_wavefront = coherent_mode_decomposition_v.get_eigenvector_wavefront(mode_index)
    output_wavefront = coherent_mode_decomposition_v.get_eigenvector_wavefront(my_mode_index)
    return output_wavefront

def run_source_gsm_h(my_mode_index=0, energy_in_keV=7):
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    if energy_in_keV == 7:
        sigma_x=3.00818e-05
        beta = 0.129748
    elif energy_in_keV == 15:
        sigma_x = 2.98958e-05
        beta = 0.0729639
    elif energy_in_keV == 30:
        sigma_x = 2.98141e-05
        beta = 0.0409362
    else:
        raise Exception("Please provide sigma_x and beta values")


    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(1e3 * energy_in_keV)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=sigma_x, amplitude=1, mode_x=0, shift=0, beta=beta)
    # previous command is useless but...
    output_wavefront.set_gaussian_hermite_mode(sigma_x=sigma_x, amplitude=1, mode_x=my_mode_index, shift=0,
                                               beta=beta)
    return output_wavefront


#
# SOURCE========================
#


def run_source_gsm_v(my_mode_index=0, energy_in_keV=7):
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    if energy_in_keV == 7:
        sigma_x=6.99408e-06
        beta=1.01172
    elif energy_in_keV == 15:
        sigma_x = 6.14502e-06
        beta = 0.624601
    elif energy_in_keV == 30:
        sigma_x = 5.73417e-06
        beta = 0.387842
    else:
        raise Exception("Please provide sigma_x and beta values")


    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-5e-05, x_max=5e-05,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(1e3 * energy_in_keV)
    # output_wavefront.set_gaussian_hermite_mode(sigma_x=5.84299e-06, amplitude=1, mode_x=0, shift=0, beta=1.56094)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=sigma_x, amplitude=1, mode_x=0, shift=0, beta=beta)
    # previous command is useless but...
    output_wavefront.set_gaussian_hermite_mode(sigma_x=sigma_x, amplitude=1, mode_x=my_mode_index, shift=0,
                                               beta=beta)
    return output_wavefront

#
# BEAMLINE========================
#


def run_beamline_h(output_wavefront,slit=50e-6, gaussian_slit=True):
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
                                       coordinates=ElementCoordinates(p=36.000000, q=0.000000,
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

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-slit/2, slit/2, -slit/2, slit/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOGaussianSlit1D, WOSlit1D
    if gaussian_slit:
        optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)
    else:
        optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    return output_wavefront


#
# BEAMLINE========================
#


def run_beamline_v(output_wavefront,slit=50e-6, gaussian_slit=True):
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
    boundary_shape = Rectangle(-slit/2, slit/2, -slit/2, slit/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOGaussianSlit1D, WOSlit1D
    if gaussian_slit:
        optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)
    else:
        optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    return output_wavefront


#
# MAIN FUNCTION========================
#


def main_h(energy_in_keV=7, source_gsm=True, slit=50e-6,gaussian_slit=False):
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(50):
        if source_gsm:
            output_wavefront = run_source_gsm_h(energy_in_keV=energy_in_keV, my_mode_index=my_mode_index)
        else:
            output_wavefront = run_source_und_h(energy_in_keV=energy_in_keV, my_mode_index=my_mode_index)
        output_wavefront = run_beamline_h(output_wavefront, slit=slit, gaussian_slit=gaussian_slit)
        tally.append(output_wavefront)

    return tally
    # tally.plot_cross_spectral_density()
    # tally.plot_spectral_density()
    # tally.plot_occupation()



def main_v(energy_in_keV=7, source_gsm=True, slit=50e-6,gaussian_slit=False):
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(50):
        if source_gsm:
            output_wavefront = run_source_gsm_v(energy_in_keV=energy_in_keV, my_mode_index=my_mode_index)
        else:
            output_wavefront = run_source_und_v(energy_in_keV=energy_in_keV, my_mode_index=my_mode_index)
        output_wavefront = run_beamline_v(output_wavefront, slit=slit, gaussian_slit=gaussian_slit)
        tally.append(output_wavefront)

    return tally
    # tally.plot_cross_spectral_density()
    # tally.plot_spectral_density()
    # tally.plot_occupation()


#
# MAIN========================
#
if __name__ == "__main__":
    import matplotlib.pylab as plt
    from srxraylib.plot.gol import plot, plot_show, set_qt
    set_qt()

    do_calculate = False
    do_plots = True

    energy_in_keV = 30 # 15 # 7



    sources = ["UND"] #, "GSM"]
    apertures = ["RECTANGULAR"] #,"GAUSSIAN"]
    subdirectory = "DataCF%d" % energy_in_keV
    #
    #
    #


    if do_calculate:
        for aperture in apertures:
            for source in sources:

                outfile = "%s/coherent_fraction_vs_slit_source_%s_aperture_%s.dat" % (subdirectory, source, aperture)

                slits = numpy.concatenate((numpy.linspace(10e-6,310e-6, 101), numpy.linspace(320e-6, 0.0015, 21)))


                f = open(outfile, "w")
                if source == "GSM":
                    source_gsm = True
                else:
                    source_gsm = False

                if aperture == "GAUSSIAN":
                    gaussian_slit = True
                else:
                    gaussian_slit = False

                for slit in slits:
                    tally_v = main_v(energy_in_keV=energy_in_keV, source_gsm=source_gsm, slit=slit, gaussian_slit=gaussian_slit)
                    tally_h = main_h(energy_in_keV=energy_in_keV, source_gsm=source_gsm, slit=slit, gaussian_slit=gaussian_slit)


                    modes_v, occ_v = tally_v.get_occupation()
                    modes_h, occ_h = tally_h.get_occupation()

                    print("slit, CF H, V: ", 1e6*slit, occ_h[0], occ_v[0])
                    f.write("%g %g %g\n" % (slit, occ_h[0], occ_v[0]))

                f.close()
                print("File written to disk: %s" % outfile)


    #
    # plots
    #
    if do_plots:
        # for source in sources:
        #     a = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_%s_aperture_GAUSSIAN.dat" % (subdirectory, source))
        #     b = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_%s_aperture_RECTANGULAR.dat" % (subdirectory, source))
        #
        #     if source == "GSM":
        #         title = "GSM source"
        #     else:
        #         title = "UNDULATOR source"
        #
        #     g = plot(
        #          1e6 * a[:,0], a[:,1],
        #          1e6 * a[:,0], a[:,2],
        #          1e6 * b[:, 0], b[:, 1],
        #          1e6 * b[:, 0], b[:, 2],
        #          legend=['Horizontal Gaussian', 'Vertical Gaussian', 'Horizontal Rectangular', 'Vertical Rectangular'],
        #          xtitle="Slit aperture [um]", ytitle="Coherent Fraction", title=title,
        #          color = ['green','blue', 'green', 'blue'],
        #          linestyle=['--','--',None,None],
        #          xlog=True, yrange=[0,1.01], show=False)
        #
        #     g[1].grid()
        #     plt.yticks(numpy.arange(0, 1.1, step=0.1))
        #
        #     plot_show()
        #
        # for aperture in apertures:
        #     b1 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_GSM_aperture_%s.dat" % (subdirectory, aperture))
        #     b2 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_UND_aperture_%s.dat" % (subdirectory, aperture))
        #
        #     g = plot(
        #         1e6 * b1[:, 0], b1[:, 1],
        #         1e6 * b1[:, 0], b1[:, 2],
        #         1e6 * b2[:, 0], b2[:, 1],
        #         1e6 * b2[:, 0], b2[:, 2],
        #         legend=['Horizontal GSM', 'Vertical GSM', 'Horizontal UND', 'Vertical UND'],
        #         xtitle="Slit aperture [um]", ytitle="Coherent Fraction", title="aperture is %s" % aperture,
        #         color=['green', 'blue', 'green', 'blue'],
        #         linestyle=['--', '--', None, None],
        #         xlog=True, yrange=[0, 1.01], show=False)
        #
        #     g[1].grid()
        #     plt.yticks(numpy.arange(0, 1.1, step=0.1))
        #
        #     plot_show()
        #
        #
        # paper
        #

        aperture = apertures[0]
        source = sources[0]
        b7_1 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_GSM_aperture_%s.dat" % ("DataCF7", aperture))
        b7_2 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_UND_aperture_%s.dat" % ("DataCF7", aperture))

        b15_1 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_GSM_aperture_%s.dat" % ("DataCF15", aperture))
        b15_2 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_UND_aperture_%s.dat" % ("DataCF15", aperture))

        b30_1 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_GSM_aperture_%s.dat" % ("DataCF30", aperture))
        b30_2 = numpy.loadtxt("%s/coherent_fraction_vs_slit_source_UND_aperture_%s.dat" % ("DataCF30", aperture))


        #
        # smooth
        #
        # from scipy.interpolate import make_interp_spline, BSpline
        # T = 1e6 * b2[:, 0]
        # power = b2[:, 1]
        # xnew = numpy.linspace(T.min(), T.max(), T.size)
        # spl = make_interp_spline(T, power, k=1)  # type: BSpline
        # power_smooth = spl(xnew)


        from scipy.signal import savgol_filter
        w7 = savgol_filter(b7_2[:, 1], 5, 2)
        w15 = savgol_filter(b15_2[:, 1], 5, 2)
        w30 = savgol_filter(b30_2[:, 1], 5, 2)



        g = plot(
            1e6 * b7_2[:, 0], w7, #b2[:, 1],
            1e6 * b7_2[:, 0], b7_2[:, 2],
            1e6 * b15_2[:, 0], w15,  # b2[:, 1],
            1e6 * b15_2[:, 0], b15_2[:, 2],
            1e6 * b30_2[:, 0], w30,  # b2[:, 1],
            1e6 * b30_2[:, 0], b30_2[:, 2],
            legend=['7 keV Horizontal', '7 keV Vertical',
                    '15 keV Horizontal', '15 keV Vertical',
                    '30 keV Horizontal', '30 keV Vertical',
                    ],
            xtitle="Slit aperture [um]", ytitle="Coherent Fraction",
            color=['green', 'green', 'blue', 'blue', 'red', 'red'],
            linestyle=[None, '--', None, '--', None, '--',],
            xlog=False, xrange=[0.,1000], yrange=[0, 1.01], show=False)

        g[1].grid()
        plt.yticks(numpy.arange(0, 1.1, step=0.1))

        def aa(x):
            return x / 565
        def invaa(x):
            return 565 * x

        # secax = g[1].secondary_xaxis('top', functions=(aa, invaa))
        # secax.set_xlabel('n')

        plt.savefig("cf_vs_aperture.eps")

        plot_show()