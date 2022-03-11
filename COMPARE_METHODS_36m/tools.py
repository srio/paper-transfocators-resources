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



def load_SRW_CSD(file, direction='x'):

    from oasys_srw.srwlib import srwl_uti_read_intens_ascii
    MI, mesh = srwl_uti_read_intens_ascii(file)

    if direction == 'x':
        CSD_Re = numpy.reshape(MI[0::2], (mesh.nx, mesh.nx))
        CSD_Im = numpy.reshape(MI[1::2], (mesh.nx, mesh.nx))
        x = numpy.linspace(mesh.xStart, mesh.xFin, mesh.nx)
        y = numpy.linspace(mesh.xStart, mesh.xFin, mesh.nx)
    else:
        CSD_Re = numpy.reshape(MI[0::2], (mesh.ny, mesh.ny))
        CSD_Im = numpy.reshape(MI[1::2], (mesh.ny, mesh.ny))
        x = numpy.linspace(mesh.yStart, mesh.yFin, mesh.ny)
        y = numpy.linspace(mesh.yStart, mesh.yFin, mesh.ny)

    CSD = CSD_Re + 1j*CSD_Im

    return CSD,x,y



def rotate_axes(csd, x1in, x2in, do_plot=0):
    x1 = x1in.copy()
    x2 = x2in.copy()
    xx1 = numpy.outer(x1, numpy.ones_like(x2))
    xx2 = numpy.outer(numpy.ones_like(x1), x2)

    X1 = numpy.linspace(x1[0], x1[-1], x1.size)
    X2 = numpy.linspace(x2[0], x2[-1], x2.size)
    XX1 = numpy.outer(X1, numpy.ones_like(X2))
    XX2 = numpy.outer(numpy.ones_like(X1), X2)

    interpolator0 = RectBivariateSpline(x1, x2, csd, bbox=[None, None, None, None], kx=3, ky=3, s=0)
    # interpolator1 = RectBivariateSpline(x1, x2, doc, bbox=[None, None, None, None], kx=3, ky=3, s=0)

    CSD = interpolator0( numpy.sqrt(2) * (XX1 + XX2) / 2, numpy.sqrt(2) * (XX2 - XX1) / 2, grid=False)
    # DOC = numpy.abs(interpolator1( numpy.sqrt(2) * (XX1 + XX2) / 2, numpy.sqrt(2) * (XX2 - XX1) / 2, grid=False))

    if do_plot:
        plot_image_with_histograms(CSD.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)sqrt(2)/2 [um]", ytitle="(x1-x2)sqrt(2)/2 [um]", title="CSD",
                                   use_profiles_instead_histograms=True)
        # plot_image_with_histograms(DOC.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)sqrt(2)/2 [um]", ytitle="(x1-x2)sqrt(2)/2 [um]", title="DoC",
        #                            use_profiles_instead_histograms=True)

    return CSD.T.copy(), X1, X2


def plotCSD(tally, rotate_axes_flag=False, srw_file=None, direction='x', range_limits=None, compare_profiles=False, normalize_to_DoC=0):

    abscissas = tally.get_abscissas()
    csd_complex = tally.get_cross_pectral_density()
    if normalize_to_DoC:
        DoC = numpy.zeros_like(csd_complex, dtype=complex)
        for i in range(abscissas.size):
            for j in range(abscissas.size):
                DoC[i,j] = csd_complex[i,j] / numpy.sqrt(numpy.abs(csd_complex[i,i]) * numpy.abs(csd_complex[j,j]))
        csd_complex = DoC

    if normalize_to_DoC:
        csd = numpy.abs(csd_complex)
    else:
        csd = numpy.abs(csd_complex)  / numpy.abs( csd_complex ).max()


    if rotate_axes_flag:
        csd, abscissas, tmp = rotate_axes(csd, abscissas, abscissas)

    if direction == 'x':
        xtitle = "x1 [um]"
        ytitle = "x2 [um]"
    else:
        xtitle = "y1 [um]"
        ytitle = "y2 [um]"

    if normalize_to_DoC:
        title="|Degree of Coherence|"
    else:
        title="|Cross Spectral Density|"
    plot_image_with_histograms(csd, abscissas * 1e6, abscissas * 1e6, title=title+" [WOFRY]",
               xtitle=xtitle, ytitle=ytitle,use_profiles_instead_histograms=True, show=0, xrange=range_limits, yrange=range_limits)


    if srw_file is not None:
        csd_srw_complex, x1, x2 = load_SRW_CSD(srw_file,direction=direction)
        if normalize_to_DoC:
            DoC_srw = numpy.zeros_like(csd_srw_complex, dtype=complex)
            for i in range(x1.size):
                for j in range(x2.size):
                    DoC_srw[i, j] = csd_srw_complex[i, j] / numpy.sqrt(
                        numpy.abs(csd_srw_complex[i, i]) * numpy.abs(csd_srw_complex[j, j]))
            csd_srw_complex = DoC_srw

        if normalize_to_DoC:
            csd_srw = numpy.abs(csd_srw_complex)
        else:
            csd_srw = numpy.abs(csd_srw_complex) / numpy.abs(csd_srw_complex).max()

        if rotate_axes_flag:
            csd_srw, x1, x2 = rotate_axes(csd_srw, x1, x2)

        plot_image_with_histograms(csd_srw, x1 * 1e6, x2 * 1e6,
                                   title=title+" [SRW]", xtitle=xtitle,
                                   ytitle=ytitle, use_profiles_instead_histograms=True, show=0,
                                   xrange=range_limits, yrange=range_limits)




    if compare_profiles:
        plot(abscissas,csd[:,csd.shape[1]//2],
             x1,csd_srw[:,csd_srw.shape[1]//2],title="H profile",legend=["WOFRY","SRW"],
                xtitle=xtitle, ytitle="", show=0)

        plot(abscissas,csd[csd.shape[0]//2,:],
             x2,csd_srw[csd_srw.shape[0]//2,:],title="V profile",legend=["WOFRY","SRW"],
                xtitle=ytitle, ytitle="", show=0)


    plot_show()


    # if new_axes:
    #     doc = csd / norm
    #     sd = numpy.sqrt(tally.get_spectral_density())
    #     norm = numpy.outer(sd,sd)
    #
    #     x1 = abscissas.copy()
    #     x2 = abscissas.copy()
    #
    #     xx1 = numpy.outer(x1, numpy.ones_like(x2))
    #     xx2 = numpy.outer(numpy.ones_like(x1), x2)
    #
    #     X1 = numpy.linspace(abscissas[0], abscissas[-1], abscissas.size)
    #     X2 = numpy.linspace(abscissas[0], abscissas[-1], abscissas.size)
    #     XX1 = numpy.outer(X1, numpy.ones_like(X2))
    #     XX2 = numpy.outer(numpy.ones_like(X1), X2)
    #
    #
    #     interpolator0 = RectBivariateSpline(x1, x2, csd, bbox=[None, None, None, None], kx=3, ky=3, s=0)
    #     interpolator1 = RectBivariateSpline(x1, x2, doc, bbox=[None, None, None, None], kx=3, ky=3, s=0)
    #
    #     CSD = numpy.abs(interpolator0((XX1 + XX2)/2, (XX2 - XX1)/2, grid=False))
    #     DOC = numpy.abs(interpolator1((XX1 + XX2)/2, (XX2 - XX1)/2, grid=False))
    #
    #     plot_image_with_histograms(CSD.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)/2 [um]", ytitle="(x1-x2)/2 [um]", title="CSD",use_profiles_instead_histograms=True)
    #     plot_image_with_histograms(DOC.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)/2 [um]", ytitle="(x1-x2)/2 [um]", title="DoC",use_profiles_instead_histograms=True)

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
        srw = numpy.loadtxt("/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/profile_I_36m_h.txt")
        srwDoC = numpy.loadtxt("/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/profile_DoC_36m_h.txt")
        fwhm_srw, quote, coordinates = get_fwhm(srw[:,1], srw[:,0])
        fwhm_srwDoC, quote, coordinates = get_fwhm(srwDoC[:, 1], srwDoC[:, 0])
        print(">>>>", srw.shape)
        plot(1e6 * abscissas, spectral_density / spectral_density.max(),
             srw[:,0],  srw[:,1] / srw[:,1].max(),
             1e6 * abscissas, y0 / y0.max(),
             numpy.sqrt(2) * srwDoC[:, 0], srwDoC[:, 1],
             legend=["Spectral Density (normalized) FWHM = %g um" % (fwhm),
                     "SRW Spectral Density (normalized) FWHM = %g um" % (fwhm_srw),
                     "Mode 0 (normalized) FWHM = %4.1f um" % (fwhm0),
                     "SRW DoC FWHM [corrected with sqrt(2)] = %g um" % (numpy.sqrt(2) * fwhm_srwDoC),],
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

def save_CSD_in_SRW_format(tally,filename="tmp.dat",direction='h'):
    abscissas = tally.get_abscissas()
    csd_complex = tally.get_cross_pectral_density()
    csd = numpy.abs( csd_complex )

    sd = numpy.sqrt(tally.get_spectral_density())
    norm = numpy.outer(sd,sd)

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
            f.write( "%g\n" % (csd_complex[i, j]).real)
            f.write( "%g\n" % (csd_complex[i, j]).imag)
    f.close()
    print("File written to disk: %s" % filename)