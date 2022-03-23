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
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import h5py
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
            number_of_points=2000,
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


def run_beamline(output_wavefront,aperture_outer=None):
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
    propagation_parameters.set_additional_parameters('magnification_x', 5.0)
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
    boundary_shape = Rectangle(-aperture_outer/2, aperture_outer/2, -aperture_outer/2, aperture_outer/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    ##########  OPTICAL ELEMENT NUMBER 3 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-(aperture_outer-5e-6)/2, (aperture_outer-5e-6)/2, -(aperture_outer-5e-6)/2, (aperture_outer-5e-6)/2)
    from wofryimpl.beamline.optical_elements.absorbers.beam_stopper import WOBeamStopper1D
    optical_element = WOBeamStopper1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    ##########  OPTICAL ELEMENT NUMBER 4 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 10 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=10.000000, q=0.000000,
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
    return output_wavefront


#
# MAIN FUNCTION========================
#


def main(aperture_outer):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(10):
        output_wavefront = run_source(my_mode_index=my_mode_index)
        output_wavefront = run_beamline(output_wavefront,aperture_outer=aperture_outer)
        tally.append(output_wavefront)

    # tally.plot_cross_spectral_density(show=1, filename="")
    # tally.plot_spectral_density(show=1, filename="")
    # tally.plot_occupation(show=1, filename="")
    return tally

from scipy.signal import savgol_filter

def get_Iav(I,ntimes=10,window=105):
    for i in range(ntimes):
        if i == 0: Iav = I.copy()
        Iav = savgol_filter(Iav, window, 3)
    return Iav

def get_fitted_DoC(I, do_plot=0):
    npoints = I.size
    Iav = get_Iav(I)
    gfit = (I / Iav) - 1
    gfinal = gfit[(npoints // 2 - npoints // 10):(npoints // 2 + npoints // 10)].max()
    print("Gfitted = ", gfinal)

    indices = numpy.arange(npoints)

    if do_plot:
        plot(
             indices, I,
             indices, Iav,
             legend=["I", "Iav"])

        plot(indices, gfit,
             indices, I * 0 + gfinal)

    print("Gfitted = ", gfinal)


    return gfinal
#
# MAIN========================
#
APERTURES = numpy.linspace(20e-6, 120e-6, 50)
CF = []
FITTED_DoC = []


for i,aperture_outer in enumerate(APERTURES):

    tally = main(aperture_outer=aperture_outer)

    x = tally.get_abscissas()
    sd = tally.get_spectral_density()
    n, occ = tally.get_occupation()
    CF.append(occ[0])
    FITTED_DoC.append(get_fitted_DoC(sd))

    # plot(x,sd, title="CF: %g aperture=%g" % (occ[0], aperture_outer*1e6) )

    if i == 0:
        SD = numpy.zeros((APERTURES.size,x.size))
    SD[i,:] = sd





plot(APERTURES*1e6,CF, xtitle="aperture(outer) [um]", ytitle="CF")
plot_image(SD, APERTURES*1e6, x*1e6, xtitle="aperture(outer) [um]",ytitle="abscissas [um]", aspect='auto')

wr = H5SimpleWriter.initialize_file(filename="wofry1d_h_doubleslit.h5", creator="srio", overwrite=1)
wr.create_entry("images", nx_default="Intensity")
wr.add_image(SD, APERTURES, x, entry_name="images", image_name="Intensity",
             title_x="distance [m]", title_y=r'X [$\mu$m]')

# f = h5py.File("wofry1d_h_doubleslit.h5", 'r')
# SD = f["images/Intensity/image_data"][()].T
# APERTURES = f["images/Intensity/axis_x"][()]
# x = f["images/Intensity/axis_y"][()]
# f.close()



filename = "wofry1d_h_doubleslit.dat"
f = open(filename,'w')
for i, aperture_outer in enumerate(APERTURES):
    f.write("%g  %g  %g  \n" % (aperture_outer, CF[i], FITTED_DoC[i]))
print("File written to disk: ", filename)