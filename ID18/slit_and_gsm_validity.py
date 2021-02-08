import numpy
from oasys.util.oasys_util import get_fwhm
from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from srxraylib.plot.gol import set_qt, plot_table, plot, plot_image
from score import Score

# see also https://people.inf.ethz.ch/gander/papers/qrneu.pdf


def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF

def run_wofry_to_screen(xmode=0, slit_size = 25e-6):
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
    plot_from_oe = 100  # set to a large number to avoid plots

    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=sigma_x, amplitude=1, mode_x=xmode, shift=0, beta=beta)

    if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

    ##########  OPTICAL SYSTEM ##########


    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle

    boundary_shape = Rectangle(-slit_size/2, slit_size/2, -slit_size/2, slit_size/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    #
    # ---- plots -----
    #
    if plot_from_oe <= 1: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 1')

    return output_wavefront


def get_CF_after_rediagonalization(sc, do_plot=False):

    # retrieve arrays
    WFs = sc.get_wavefronts()
    nmodes = sc.get_number_of_calls()
    abscissas = WFs[-1].get_abscissas()

    #
    # calculate the CSD
    #

    input_array = numpy.zeros((nmodes, abscissas.size), dtype=complex)
    for i,wf in enumerate(WFs):
        input_array[i,:] = wf.get_complex_amplitude() # tmp[i][0]

    cross_spectral_density = numpy.zeros((abscissas.size, abscissas.size), dtype=complex)

    for i in range(nmodes):
        cross_spectral_density += numpy.outer(numpy.conjugate(input_array[i, :]), input_array[i, :])

    if do_plot:
        plot_image(numpy.abs(cross_spectral_density))
    print("matrix cross_spectral_density: ", cross_spectral_density.shape)

    #
    # diagonalize the CSD
    #

    w, v = numpy.linalg.eig(cross_spectral_density)
    print(w.shape, v.shape)
    idx = w.argsort()[::-1]  # large to small
    eigenvalues  = numpy.real(w[idx])
    eigenvectors = v[:, idx].T



    #
    # plot intensity
    #
    if do_plot:
        y = numpy.zeros_like(abscissas)

        for i in range(nmodes):
            y += eigenvalues[i] * numpy.real(numpy.conjugate(eigenvectors[i,:]) * eigenvectors[i,:])

        cumulated_intensity = sc.get_additional_stored_values()[-1][0]

        plot(abscissas, cumulated_intensity,
             abscissas, y, legend=["Data", "From modes"])

        print(nmodes, eigenvalues.shape)

        plot(numpy.arange(nmodes), eigenvalues[0:nmodes]/(eigenvalues[0:nmodes].sum()),
             title="mode occupation")

    return eigenvalues[0] / eigenvalues.sum()

def get_gsm_cf(factor=1.0):
    gsm_slit = slit_size / factor
    gsm_sigma_x = numpy.sqrt(1.0 / (1.0 / sigma_x ** 2 + 1.0 / gsm_slit ** 2))
    gsm_xi_initial = sigma_x * beta
    gsm_xi = numpy.sqrt(1.0 / (1.0 / gsm_xi_initial ** 2 + 1.0 / gsm_slit ** 2))
    return get_coherent_fraction_exact(gsm_xi / gsm_sigma_x)

if __name__ == "__main__":

    #
    # inputs
    #
    beta = 1.151 #0.0922395 #0.02 #
    sigma_x = 3.03783e-05
    slit_size_over_sigma_x = 10.0

    SLIT_SIZE_OVER_SIGMA = numpy.linspace(0.5,10,30)
    CF = numpy.zeros_like(SLIT_SIZE_OVER_SIGMA)
    GSM_CF_4 = numpy.zeros_like(SLIT_SIZE_OVER_SIGMA)
    GSM_CF_12 = numpy.zeros_like(SLIT_SIZE_OVER_SIGMA)
    GSM_CF_235 = numpy.zeros_like(SLIT_SIZE_OVER_SIGMA)
    #
    #
    #
    set_qt()

    for i, slit_size_over_sigma_x in enumerate(SLIT_SIZE_OVER_SIGMA):
        slit_size = slit_size_over_sigma_x * sigma_x

        sc = Score(scan_variable_name='mode index',
                   additional_stored_variable_names=['cumulated_intensity'],
                   do_store_wavefronts=True)

        for xmode in range(51):

            # output_wavefront = run_wofry_source(xmode=xmode)
            output_wavefront = run_wofry_to_screen(xmode=xmode, slit_size=slit_size)

            if xmode == 0:
                intens = output_wavefront.get_intensity()
            else:
                intens += output_wavefront.get_intensity()

            sc.append(output_wavefront,
                      scan_variable_value=xmode,
                      additional_stored_values=[intens])


        cf = get_CF_after_rediagonalization(sc, do_plot=0)
        cf_source = get_coherent_fraction_exact(beta)

        CF[i] = cf
        # factor = numpy.sqrt(12)
        # factor = 2.355
        # factor = 4.1


        # gsm_slit = slit_size /  factor
        # gsm_sigma_x = numpy.sqrt(1.0 / (1.0 / sigma_x**2 + 1.0 / gsm_slit**2))
        # gsm_xi_initial = sigma_x * beta
        # gsm_xi = numpy.sqrt(1.0 / (1.0 / gsm_xi_initial**2 + 1.0 / gsm_slit**2))
        # gsm_cf = get_coherent_fraction_exact(gsm_xi / gsm_sigma_x)

        GSM_CF_12[i] = get_gsm_cf(factor=numpy.sqrt(12))
        GSM_CF_4[i] = get_gsm_cf(factor=4.1)
        GSM_CF_235[i] = get_gsm_cf(factor=2.355)


        print("Exact CF at the source: ", cf_source)
        print("Slit size: %g um" % (slit_size * 1e6))
        print("CF (rediag): %g " % cf)
        print("CF/CFsource: %g " % (cf / cf_source))
        # print("CF (gsm): %g " % gsm_cf)
        # print("CF/CFsource: %g " % (gsm_cf / cf_source))

        # CF[i] = cf
        # GSM_CF[i] = gsm_cf

    plot(SLIT_SIZE_OVER_SIGMA, CF,
         SLIT_SIZE_OVER_SIGMA, GSM_CF_4,
         SLIT_SIZE_OVER_SIGMA, GSM_CF_12,
         SLIT_SIZE_OVER_SIGMA, GSM_CF_235,
         xtitle="slit size over sigma", ytitle="Coherent Fraction",
         title="beta: %g, sigma_x: %g um" % (beta, sigma_x*1e6),
         legend=["Exact",
                 "Approx GSM W/4.1",
                 "Approx GSM W/sqrt(12)",
                 "Approx GSM W/2.355",])



