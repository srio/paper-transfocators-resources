from srxraylib.plot.gol import plot, plot_image
import numpy

def run_wofry(scan_direction='H'):
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


    plot_from_oe = 0 # set to a large number to avoid plots


    ##########  SOURCE ##########


    #
    # create output_wavefront
    #
    #
    if scan_direction == 'H':
        sigmaxx=2.97321e-05
        sigmaxpxp=4.37237e-06
    else:
        sigmaxx = 5.2915e-06
        sigmaxpxp = 1.88982e-06

    from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import UndulatorCoherentModeDecomposition1D
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
        scan_direction=scan_direction,
        sigmaxx=sigmaxx,
        sigmaxpxp=sigmaxpxp,
        useGSMapproximation=False,)
    # make calculation
    coherent_mode_decomposition_results = coherent_mode_decomposition.calculate()

    mode_index = 0
    output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(mode_index)

    return coherent_mode_decomposition_results, output_wavefront

if __name__ == "__main__":
    coherent_mode_decomposition_results, output_wavefront  =run_wofry(scan_direction='V')

    # if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')


    csd = numpy.abs(coherent_mode_decomposition_results["CSD"])
    abscissas = coherent_mode_decomposition_results["abscissas"]
    plot_image(csd, abscissas * 1e6, abscissas * 1e6, title="Cross spectral density", xtitle="x1 [um]",
               ytitle="x2 [um]")

    x1 = abscissas.copy()
    x2 = abscissas.copy()

    xx1 = numpy.outer(x1, numpy.ones_like(x2))
    xx2 = numpy.outer(numpy.ones_like(x1), x2)

    X1 = numpy.linspace(abscissas[0], abscissas[-1], abscissas.size)
    X2 = numpy.linspace(abscissas[0], abscissas[-1], abscissas.size)
    XX1 = numpy.outer(X1, numpy.ones_like(X2))
    XX2 = numpy.outer(numpy.ones_like(X1), X2)

    from scipy.interpolate import interp2d
    f = interp2d(x1, x2, numpy.abs(csd), kind='linear')

    # CSD = numpy.zeros((X1.size,X2.size))
    #
    # for i in range(X1.size):
    #     for j in range(X2.size):
    #         CSD[i,j] = f((X1[i]+X2[j]),(X2[j]-X1[i]))


    from scipy.interpolate import RectBivariateSpline
    interpolator0 = RectBivariateSpline(x1, x2, csd, bbox=[None, None, None, None], kx=3, ky=3, s=0)

    CSD = numpy.abs(interpolator0((XX1 + XX2), (XX2 - XX1), grid=False)) #/ \
                    # numpy.sqrt(numpy.abs(interpolator0(X + Y, X + Y, grid=False))) / \
                    # numpy.sqrt(numpy.abs(interpolator0(XX2 - XX1, XX2 - XX1, grid=False)))


    # CSD = f(X1, X2)
    #
    plot_image(CSD, X1 * 1e6, X2 * 1e6, title="", xtitle="x1-x2 [um]", ytitle="x1+x2 [um]")

    profile = CSD[:,X2.size//2]
    mode0 = numpy.abs(output_wavefront.get_complex_amplitude()) ** 2
    indices = numpy.arange(x1.size)
    intensity = csd[indices,indices]
    profileI = CSD[X2.size//2, :]
    profile_vs_x1 = csd[: , x2.size//2]

    plot(1e6 * X1, profile / profile.max(),
         1e6 * X1, profileI / profileI.max(),
         1e6 * abscissas, (mode0 / mode0.max()),
         1e6 * abscissas, (intensity / intensity.max()),
         1e6 * abscissas, (profile_vs_x1 / profile_vs_x1.max()),
         legend=['I vs (x2-x1)/2','I vs (x1+x2)/2','mode1','intensity','csd(x1,0)'], xrange=[-100,100])

    # interpolator0 = RectBivariateSpline(coor, coor_conj, mutual_intensity, bbox=[None, None, None, None], kx=3, ky=3, s=0)
    #
    # X = numpy.outer(coor, numpy.ones_like(coor_conj))
    # Y = numpy.outer(numpy.ones_like(coor), coor_conj)
    #
    # nmResDegCoh_z = numpy.abs(interpolator0(X + Y, X - Y, grid=False)) / \
    #                 numpy.sqrt(numpy.abs(interpolator0(X + Y, X + Y, grid=False))) / \
    #                 numpy.sqrt(numpy.abs(interpolator0(X - Y, X - Y, grid=False)))