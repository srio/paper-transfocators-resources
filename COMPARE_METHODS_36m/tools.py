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

import pylab as plt


def load_SRW_SD(file_name):

    from oasys_srw.srwlib import srwl_uti_read_intens_ascii
    MI, mesh = srwl_uti_read_intens_ascii(file_name)

    SD = numpy.reshape(MI,(mesh.ny, mesh.nx)).T

    x = numpy.linspace(mesh.xStart,mesh.xFin, mesh.nx)
    y = numpy.linspace(mesh.yStart,mesh.yFin, mesh.ny)

    return SD, x, y



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



def rotate_axes(csd, x1in, x2in,
                rotate_axes_normalization = 0,  # 0=sum, 1=half-sum [SRW], 2=rotated
                do_plot=0):
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


    if rotate_axes_normalization == 0:
        CSD = interpolator0((XX1 + XX2), (XX2 - XX1), grid=False)
    elif rotate_axes_normalization == 1:
        CSD = interpolator0(numpy.sqrt(1) * (XX1 + XX2) / 2, numpy.sqrt(1) * (XX2 - XX1) / 2, grid=False)
    elif rotate_axes_normalization == 2:
        CSD = interpolator0(numpy.sqrt(2) * (XX1 + XX2) / 2, numpy.sqrt(2) * (XX2 - XX1) / 2, grid=False)
    elif rotate_axes_normalization == 3:
        CSD = interpolator0( (XX1 + XX2) / 2, (XX2 - XX1) , grid=False)

    # DOC = numpy.abs(interpolator1( numpy.sqrt(2) * (XX1 + XX2) / 2, numpy.sqrt(2) * (XX2 - XX1) / 2, grid=False))

    if do_plot:
        plot_image_with_histograms(CSD.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)sqrt(2)/2 [um]", ytitle="(x1-x2)sqrt(2)/2 [um]", title="CSD",
                                   use_profiles_instead_histograms=True)
        # plot_image_with_histograms(DOC.T, X2 * 1e6, X1 * 1e6, xtitle="(x1+x2)sqrt(2)/2 [um]", ytitle="(x1-x2)sqrt(2)/2 [um]", title="DoC",
        #                            use_profiles_instead_histograms=True)

    return CSD.T.copy(), X1, X2

def rotate_axes_new(csd, x1in, x2in,
                rotate_axes_normalization = 0,  # 0=sum, diff, 1=half-sum, half-diff [SRW], 2=rotated, 3=half-sum, diff
                do_plot=0):

    x1 = x1in.copy()
    x2 = x2in.copy()
    interpolator0 = RectBivariateSpline(x1, x2, csd, bbox=[None, None, None, None], kx=3, ky=3, s=0)


    if rotate_axes_normalization == 0:
        factor_x = factor_y = 1.0
    elif rotate_axes_normalization == 1:
        factor_x = factor_y = 0.5
    elif rotate_axes_normalization == 2:
        factor_x = factor_y = 1.0 / numpy.sqrt(2)
    elif rotate_axes_normalization == 3:
        factor_x = 0.5
        factor_y = 1.0


    x1m = x1.min() #  x1[0]
    x1M = x1.max()  #  x1[-1]
    x2m = x2.min() #  x2[0]
    x2M = x2.max()  #  x2[-1]

    print(">>>> neww ORIGINAL extrema: ", x1m, x1M, x2m, x2M)
    X1 = numpy.linspace(factor_x * (x1m + x2m), factor_x * (x1M + x2M), x1.size)
    X2 = numpy.linspace(factor_y * (-x1M + x2m), factor_y * (x2M - x1m), x2.size)

    print(">>>> neww X1 extrema: ", X1[0], X1[-1], X1.min(), X1.max())
    print(">>>> neww X2 extrema: ", X2[0], X2[-1], X2.min(), X2.max())

    print(">>>> neww limits H: ", factor_x * (x1m + x2m), factor_y * (x1M + x2M))
    print(">>>> neww limits V: ", factor_y * (x1m - x2M), factor_y * (x1M - x2m))
    XX1 = numpy.outer(X1, numpy.ones_like(X2))
    XX2 = numpy.outer(numpy.ones_like(X1), X2)

    # interpolator1 = RectBivariateSpline(x1, x2, doc, bbox=[None, None, None, None], kx=3, ky=3, s=0)


    print(">>>> shapes", XX1.shape, XX2.shape)
    CSD = interpolator0( 0.5 * (XX1/factor_x - X2/factor_y), 0.5 * (XX2/factor_y + XX1/factor_x), grid=False)
    print(">>>> shapes", CSD.shape)


    return CSD.copy(), X1, X2




def plotCSD(tally,
            rotate_axes_flag=0,  # 0=No, 1=OLD!!!!!!!!!! NOT GOOD!,  2=New
            rotate_axes_normalization=0, # 0=sum, diff, 1=half-sum, half-diff [SRW], 2=rotated, 3=half-sum, diff
            srw_file=None,
            direction='x',
            range_limits=None,
            compare_profiles=0, # 0=No, 1=Yes, 2=Yes plus first mode
            normalize_to_DoC=0,
            do_plot=0):

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

    print(">>>> shapes BEFORE", csd.shape, abscissas.shape)
    if rotate_axes_flag == 0:
        abscissas1 = abscissas
        abscissas2 = abscissas
    elif rotate_axes_flag == 1:
        csd, abscissas1, abscissas2 = rotate_axes(csd, abscissas, abscissas, rotate_axes_normalization=rotate_axes_normalization)
    elif rotate_axes_flag == 2:
        csd, abscissas1, abscissas2 = rotate_axes_new(csd, abscissas, abscissas,
                                                      rotate_axes_normalization=rotate_axes_normalization)

    print(">>>> shapes AFTER", csd.shape, abscissas2.shape)
    if direction == 'x':
        if rotate_axes_flag:
            if rotate_axes_normalization == 0:
                xtitle = "$(x_1+x_2)$ [$\mu$m]"
                ytitle = "$(x_2-x_1)$ [$\mu$m]"
            elif rotate_axes_normalization == 1:
                xtitle = "$(x_1+x_2)/2$ [$\mu$m]"
                ytitle = "$(x_2-x_1)/2$ [$\mu$m]"
            elif rotate_axes_normalization == 2:
                xtitle = "$(x_1+x_2)/\sqrt{2}$ [$\mu$m]"
                ytitle = "$(x_2-x_1)/\sqrt{2}$ [$\mu$m]"
            elif rotate_axes_normalization == 3:
                xtitle = "$(x_1+x_2)/2$ [$\mu$m]"
                ytitle = "$(x_2-x_1)$ [$\mu$m]"

        else:
            xtitle = "$x_1$ [$\mu$m]"
            ytitle = "$x_2$ [$\mu$m]"

        xtitleProfile = "$x_1$ [$\mu$m]"
        ytitleProfile = "$x_2$ [$\mu$m]"
    else:
        if rotate_axes_flag:
            if rotate_axes_normalization == 0:
                xtitle = "$(x_1+x_2)$ [$\mu$m]"
                ytitle = "$(x_2-x_1)$ [$\mu$m]"
            elif rotate_axes_normalization == 1:
                xtitle = "$(x_1+x_2)/2$ [$\mu$m]"
                ytitle = "$(x_2-x_1)/2$ [$\mu$m]"
            elif rotate_axes_normalization == 2:
                xtitle = "$(x_1+x_2)/\sqrt{2}$ [$\mu$m]"
                ytitle = "$(x_2-x_1)/\sqrt{2}$ [$\mu$m]"
            elif rotate_axes_normalization == 3:
                xtitle = "$(x_1+x_2)/2$ [$\mu$m]"
                ytitle = "$(x_2-x_1)$ [$\mu$m]"
        else:
            xtitle = "$y_1$ [$\mu$m]"
            ytitle = "$y_2$ [$\mu$m]"

        xtitleProfile = "$y_1$ [$\mu$m]"
        ytitleProfile = "$y_2$ [$\mu$m]"

    if normalize_to_DoC:
        title="|Degree of Coherence|"
    else:
        title="|Cross Spectral Density|"

    if do_plot == 1:
        plot_image(csd, abscissas1 * 1e6, abscissas2 * 1e6, title=title+" [WOFRY]",
                   xtitle=xtitle, ytitle=ytitle, show=0, xrange=range_limits, yrange=range_limits)
    elif do_plot == 2:
        plot_image_with_histograms(csd, abscissas1 * 1e6, abscissas2 * 1e6, title=title+" [WOFRY]",
                   xtitle=xtitle, ytitle=ytitle, show=0, xrange=range_limits, yrange=range_limits,use_profiles_instead_histograms=True,)

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

        if rotate_axes_flag == 0:
            pass
        elif rotate_axes_flag == 1:
            csd_srw, x1, x2 = rotate_axes(csd_srw, x1, x2, rotate_axes_normalization=rotate_axes_normalization)
        elif rotate_axes_flag == 2:
            csd_srw, x1, x2 = rotate_axes_new(csd_srw, x1, x2, rotate_axes_normalization=rotate_axes_normalization)



        if do_plot == 1:
            plot_image(csd_srw, x1 * 1e6, x2 * 1e6,
                                       title=title+" [SRW]", xtitle=xtitle,
                                       ytitle=ytitle, show=0,
                                       xrange=range_limits, yrange=range_limits)
        elif do_plot == 2:
            plot_image_with_histograms(csd_srw, x1 * 1e6, x2 * 1e6,
                                       title=title+" [SRW]", xtitle=xtitle,
                                       ytitle=ytitle, show=0,
                                       xrange=range_limits, yrange=range_limits, use_profiles_instead_histograms=True,)

    if srw_file is not None and compare_profiles > 0:
        if do_plot > 0:
            plot(1e6*abscissas1, csd[:, csd.shape[1] // 2],
                 1e6*x1, csd_srw[:, csd_srw.shape[1] // 2], title="H profile", legend=[
                    "WOFRY (fwhm: %g) " % (fwhm(1e6*abscissas1, csd[:, csd.shape[1] // 2], xrange=range_limits)),
                    "SRW (fwhm: %g) "   % (fwhm(1e6*x1, csd_srw[:, csd_srw.shape[1] // 2], xrange=range_limits)),
                ],
                 xtitle=xtitleProfile, ytitle="", xrange=range_limits, yrange=[0,1.1], show=0)

        if compare_profiles == 2:
            eigenvalues = tally.get_eigenvalues()
            eigenvectors = tally.get_eigenvectors()
            y0 = eigenvalues[0] * numpy.real(numpy.conjugate(eigenvectors[0, :]) * eigenvectors[0, :])

            if do_plot > 0:
                plot(
                     1e6*abscissas2, csd[csd.shape[0] // 2, :],
                     1e6*x2, csd_srw[csd_srw.shape[0] // 2, :],
                     1e6*abscissas, y0/y0.max(),
                     title="V profile", legend=[
                        "WOFRY (fwhm: %g) "  % (fwhm(1e6*abscissas2, csd[csd.shape[0] // 2, :], xrange=range_limits)),
                        "SRW (fwhm: %g) "    % (fwhm(1e6*x2, csd_srw[csd_srw.shape[0] // 2, :], xrange=range_limits)),
                        "Mode 0 (fwhm: %g) " % (fwhm(1e6*abscissas, y0/y0.max())),
                    ],
                     xtitle=ytitleProfile, ytitle="", xrange=range_limits, yrange=[0,1.1], show=0)

        elif compare_profiles == 1:
            if do_plot:
                plot(1e6*abscissas2,csd[csd.shape[0]//2,:],
                     1e6*x2,csd_srw[csd_srw.shape[0]//2,:],title="V profile",legend=[
                        "WOFRY (fwhm: %g)" % (fwhm(1e6*abscissas2,csd[csd.shape[0]//2,:], xrange=range_limits)),
                        "SRW (fwhm: %g)"   % (fwhm(1e6*x2,csd_srw[csd_srw.shape[0]//2,:], xrange=range_limits)),
                    ],
                        xtitle=ytitleProfile, ytitle="", xrange=range_limits, yrange=[0,1.1], show=0)


    if do_plot:
        plot_show()

    return [csd, csd_srw],[abscissas1*1e6, x1*1e6], [abscissas2*1e6, x2*1e6]


def plot_four_images_with_histograms(z_list,x_list,y_list,
                               title_list=None, xtitle=r"X", ytitle=r"Y",
                               xrange=None, yrange=None,
                               cmap=None, aspect='auto', show=True,
                               add_colorbar=False, figsize=(8,8),
                               use_profiles_instead_histograms=False,
                               ):

    if aspect is None: aspect == 'auto'

    figure = plt.figure(figsize=figsize)

    SHIFT_H = [0,    0.5, 0, 0.5]
    SHIFT_V = [0.5,  0.5, 0, 0]

    for i in range(len(z_list)):
        z = z_list[i].T
        x = x_list[i]
        y = y_list[i]
        shift_h = SHIFT_H[i]
        shift_v = SHIFT_V[i]



        if xrange is None:
            xrange = [x.min(),x.max()]

        if yrange is None:
            yrange = [y.min(),y.max()]


        hfactor = 1.0
        vfactor = 1.0

        # left, width = 0.1, 0.6
        # bottom, height = 0.1, 0.6
        # bottom_h = left_h = left + width + 0.02

        # rect_scatter = [left, bottom, width, height]
        # rect_histx   = [left, bottom_h, width, 0.2]
        # rect_histy   = [left_h, bottom, 0.2, height]



        width  = 0.6 / 2
        height = 0.6 / 2

        left   = 0.1 + shift_h
        bottom = 0.1 + shift_v

        left_h = left + width + 0.02
        bottom_h = bottom + height + 0.02
        hh = 0.2 / 4

        rect_scatter = [left, bottom, width, height]
        rect_histx   = [left, bottom_h, width, hh]
        rect_histy   = [left_h, bottom, hh, height]

        #
        #main plot
        #
        axScatter = figure.add_axes(rect_scatter)

        if isinstance(xtitle, list):
            axScatter.set_xlabel(xtitle[i])
        else:
            axScatter.set_xlabel(xtitle)

        if isinstance(ytitle, list):
            axScatter.set_ylabel(ytitle[i])
        else:
            axScatter.set_ylabel(ytitle)

        if aspect == 'equal':
            axScatter.set_aspect(aspect)

        axScatter.axis(xmin=hfactor*xrange[0],xmax=xrange[1])
        axScatter.axis(ymin=vfactor*yrange[0],ymax=yrange[1])



        axs = axScatter.pcolormesh(x,y,z,cmap=cmap)

        #
        #histograms
        #

        if aspect == 'equal':
            pos0 = axScatter.get_position()
            mm = np.min((pos0.height, pos0.width)) * 0.6
            axHistx = figure.add_axes([pos0.x0, pos0.y0 +pos0.height, pos0.width, mm], sharex=axScatter)
            axHisty = figure.add_axes([pos0.x0 + pos0.width, pos0.y0, mm * figsize[1] / figsize[0], pos0.height], sharey=axScatter)
        else:
            axHistx = figure.add_axes(rect_histx, sharex=axScatter)
            axHisty = figure.add_axes(rect_histy, sharey=axScatter)

        if use_profiles_instead_histograms:
            hx = z[z.shape[0]//2, :]
            hy = z[:, z.shape[1]//2]
        else:
            hx = z.sum(axis=0)
            hy = z.sum(axis=1)

        axHistx.plot(x,hx)
        axHisty.plot(hy,y)

        # supress ordinates labels ans ticks
        axHistx.get_yaxis().set_visible(False)
        axHisty.get_xaxis().set_visible(False)
        axHistx.set_ylim(ymin=0, ymax=hy.max()*1.1)
        axHisty.set_xlim(xmin=0, xmax=hx.max() * 1.1)
        # supress abscissas labels (keep ticks)
        for tl in axHistx.get_xticklabels(): tl.set_visible(False)
        for tl in axHisty.get_yticklabels(): tl.set_visible(False)

        if title_list is not None:
            axHistx.set_title(title_list[i])


    if add_colorbar:
        plt.colorbar(axs)

    if show:
        plt.show()

    return figure, axScatter, axHistx, axHisty

#
# MAIN FUNCTION========================
#

def fwhm(x,y, xrange=None):
    mask = x * 0 + 1.0
    if xrange is not None:
        for i in range(x.size):
            if x[i] < xrange[0]:
                mask[i] = 0.0
            if x[i] > xrange[1]:
                mask[i] = 0.0
    fwhm1, quote, coordinates = get_fwhm(y * mask, x)
    return fwhm1

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