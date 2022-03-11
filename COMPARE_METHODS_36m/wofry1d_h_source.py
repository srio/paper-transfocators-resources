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

from oasys.util.oasys_util import get_fwhm
from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show
from scipy.interpolate import RectBivariateSpline

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
            number_of_points=1000,
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


def run_beamline(output_wavefront):
    return output_wavefront

def main(do_plot=1):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(50):
        output_wavefront = run_source(my_mode_index=my_mode_index)
        output_wavefront = run_beamline(output_wavefront)
        tally.append(output_wavefront)

    # tally.plot_cross_spectral_density(show=1, filename="")
    # tally.plot_spectral_density(show=1, filename="")
    # tally.plot_occupation(show=1, filename="")

    return tally


#
# MAIN========================
#

from tools import plotCSD

tally = main(do_plot=0)
plotCSD(
    tally, range_limits=[-100,100], compare_profiles=2, rotate_axes_flag=1, normalize_to_DoC=1,
    direction='x',
    srw_file='/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/source_imaging/CSDx7.0keV_50k_ME_intensity.dat',
    # srw_file='/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/source_imaging/CSDx7.0keV_15k_ME_intensity.dat',
      )