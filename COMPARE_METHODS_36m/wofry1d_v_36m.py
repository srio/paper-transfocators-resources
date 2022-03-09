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
from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms
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
            number_of_points=800,
            distance_to_screen=100,
            scan_direction='V',
            sigmaxx=5.2915e-06,
            sigmaxpxp=1.88982e-06,
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
    return output_wavefront


def plot2(tally):

    abscissas = tally.get_abscissas()
    csd = numpy.abs( tally.get_cross_pectral_density() )

    sd = numpy.sqrt(tally.get_spectral_density())
    norm = numpy.outer(sd,sd)
    doc = csd / norm



    plot_image_with_histograms(csd, abscissas * 1e6, abscissas * 1e6, title="Cross spectral density", xtitle="x1 [um]",
               ytitle="x2 [um]",use_profiles_instead_histograms=True)

    x1 = abscissas.copy()
    x2 = abscissas.copy()

    xx1 = numpy.outer(x1, numpy.ones_like(x2))
    xx2 = numpy.outer(numpy.ones_like(x1), x2)

    X1 = numpy.linspace(abscissas[0], abscissas[-1], abscissas.size)
    X2 = numpy.linspace(abscissas[0], abscissas[-1], abscissas.size)
    XX1 = numpy.outer(X1, numpy.ones_like(X2))
    XX2 = numpy.outer(numpy.ones_like(X1), X2)


    interpolator0 = RectBivariateSpline(x1, x2, csd, bbox=[None, None, None, None], kx=3, ky=3, s=0)
    interpolator1 = RectBivariateSpline(x1, x2, doc, bbox=[None, None, None, None], kx=3, ky=3, s=0)

    CSD = numpy.abs(interpolator0((XX1 + XX2)/2, (XX2 - XX1)/2, grid=False))
    DOC = numpy.abs(interpolator1((XX1 + XX2)/2, (XX2 - XX1)/2, grid=False))

    plot_image_with_histograms(CSD.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)/2 [um]", ytitle="(x1-x2)/2 [um]", title="CSD",use_profiles_instead_histograms=True)
    plot_image_with_histograms(DOC.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)/2 [um]", ytitle="(x1-x2)/2 [um]", title="DoC",use_profiles_instead_histograms=True)

    # profile = CSD[:,X2.size//2]
    # mode0 = numpy.abs(output_wavefront.get_complex_amplitude()) ** 2
    # indices = numpy.arange(x1.size)
    # intensity = csd[indices,indices]
    # profileI = CSD[X2.size//2, :]
    # profile_vs_x1 = csd[: , x2.size//2]
    #
    # plot(1e6 * X1, profile / profile.max(),
    #      1e6 * X1, profileI / profileI.max(),
    #      1e6 * abscissas, (mode0 / mode0.max()),
    #      1e6 * abscissas, (intensity / intensity.max()),
    #      1e6 * abscissas, (profile_vs_x1 / profile_vs_x1.max()),
    #      legend=['I vs (x2-x1)/2','I vs (x1+x2)/2','mode1','intensity','csd(x1,0)'], xrange=[-100,100])


#
# MAIN FUNCTION========================
#

def plot1(tally, add_srw=0):

    abscissas = tally.get_abscissas()
    eigenvalues = tally.get_eigenvalues()
    eigenvectors = tally.get_eigenvectors()

    spectral_density = tally.get_spectral_density() # numpy.zeros_like(abscissas)
    fwhm, quote, coordinates = get_fwhm(spectral_density, 1e6 * abscissas)

    y0 = eigenvalues[0] * numpy.real(numpy.conjugate(eigenvectors[0, :]) * eigenvectors[0, :])
    fwhm0, quote, coordinates = get_fwhm(y0, 1e6 * abscissas)

    if add_srw:
        srw = numpy.loadtxt("/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/profile_I_36m_v.txt")
        srwDoC = numpy.loadtxt("/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/profile_DoC_36m_v.txt")
        fwhm_srw, quote, coordinates = get_fwhm(srw[:,1], srw[:,0])
        fwhm_srwDoC, quote, coordinates = get_fwhm(srwDoC[:, 1], srwDoC[:, 0])
        plot(1e6 * abscissas, spectral_density / spectral_density.max(),
             srw[:,0],  srw[:,1] / srw[:,1].max(),
             1e6 * abscissas, y0 / y0.max(),
             numpy.sqrt(2) * srwDoC[:, 0], srwDoC[:, 1] ,
             legend=["Spectral Density (normalized) FWHM = %g um" % (fwhm),
                     "SRW Spectral Density (normalized) FWHM = %g um" % (fwhm_srw),
                     "Mode 0 (normalized) FWHM = %4.1f um" % (fwhm0),
                     "SRW DoC FWHM [corrected with sqrt(2)] =  %4.1f um" % (numpy.sqrt(2) * fwhm_srwDoC),],
             xtitle="x [um]", ytitle="(a.u)", show=True, xrange=[-1000,1000])
    else:
        plot(1e6 * abscissas, spectral_density / spectral_density.max(),
             1e6 * abscissas, y0 / y0.max(),
             legend=["Spectral Density (normalized) FWHM = %g um" % (fwhm),
                     "Mode 0 (normalized) FWHM = %4.1f um" % (fwhm0), ],
             xtitle="x [um]", ytitle="(a.u)", show=True)


    # csd = tally.get_cross_pectral_density()
    # plot_image(numpy.abs(csd), 1e6 * tally.abscissas, 1e6 * tally.abscissas,
    #            title="Cross Spectral Density", xtitle="X1 [um]", ytitle="X2 [um]", show=True)

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

def save_CSD_in_SRW_format(tally,filename="tmp.dat",direction='h'):
    abscissas = tally.get_abscissas()
    csd_complex = tally.get_cross_pectral_density()
    csd = numpy.abs( csd_complex )

    sd = numpy.sqrt(tally.get_spectral_density())
    norm = numpy.outer(sd,sd)
    doc = csd / norm



    plot_image_with_histograms(csd, abscissas * 1e6, abscissas * 1e6, title="Cross spectral density", xtitle="x1 [um]",
               ytitle="x2 [um]",use_profiles_instead_histograms=True)


    f = open(filename, 'w')

    f.write("# Complex Mutual Intensity [ph/s/.1%bw/mm^2] (C-aligned, inner loop is vs Photon Energy, outer loop vs Vertical Position)\n")
    f.write("# 7000.0 #Initial Photon Energy [eV]\n")
    f.write("# 7000.0 #Final Photon Energy [eV]\n")
    f.write("# 1 #Number of points vs Photon Energy\n")
    if direction == 'h':
        f.write("# %g #Initial Horizontal Position [m]\n" % abscissas[0])
        f.write("# %g #Final Horizontal Position [m]\n" % abscissas[-1])
        f.write("# %d #Number of points vs Horizontal Position\n" % (abscissas.size))
        f.write("# 0 #Initial Vertical Position [m]\n")
        f.write("# 0 #Final Vertical Position [m]\n")
        f.write("# 1 #Number of points vs Vertical Position\n")
    else:
        f.write("# 0 #Initial Horizontal Position [m]\n")
        f.write("# 0 #Final Horizontal Position [m]\n")
        f.write("# 1 #Number of points vs Horizontal Position\n")
        f.write("# %g #Initial Vertical Position [m]\n" % (abscissas[0]))
        f.write("# %g #Final Vertical Position [m]\n" % (abscissas[-1]))
        f.write("# %d #Number of points vs Vertical Position\n" % (abscissas.size))
    f.write("# 1 #Number of components\n")

    for j in range(abscissas.size):
        for i in range(abscissas.size):
            f.write(repr(csd_complex[i,j])+"\n")
    f.close()
    print("File written to disk: %s" % filename)
#
# MAIN========================
#


tally = main(do_plot=0)
save_CSD_in_SRW_format(tally, filename="csd_v.dat", direction='v')
plot1(tally, add_srw=1)
plot2(tally)