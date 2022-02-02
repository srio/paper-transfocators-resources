#!/usr/bin/python
# coding: utf-8

import barc4plots.barc4plots as b4pt
import logging.handlers
import datetime
import time
import os

from numpy import reshape, linspace, zeros, trapz, amax, argmax, interp, absolute
from scipy.interpolate import interp2d
from scipy.ndimage import gaussian_filter

import sys
sys.path.insert(0, '../srw_python')

from srwlib import srwl_uti_save_intens_ascii, srwl_uti_read_intens_ascii, SRWLStokes


def get_slice(_image, _x, _y, _coords_x, _coords_y, _integrate=False):
    if _integrate is False:
        f = interp2d(_x, _y, _image, kind='linear')

        if _coords_x == ':':
            # coords_x = x
            _integrated = zeros(_x.size)
            _cut = f(_x, _coords_y)
            _axis = _x
            # return _cut, _x
        if _coords_y == ':':
            # _coords_y = y
            _integrated = zeros(_y.size)
            _cut = f(_coords_x, _y)
            _axis = _y

        for i in range(_integrated.size):
            _integrated[i] = _cut[i]
        return _integrated, _axis
    else:
        ny = _y.size
        nx = _x.size
        if _coords_x == ':':
            _integrated = zeros(nx)
            for i in range(nx):
                _integrated[i] = trapz(_image[:,i])
            return _integrated, _x
        if _coords_y == ':':
            _integrated = zeros(ny)
            for i in range(ny):
                _integrated[i] = trapz(_image[i,:])
            return _integrated, _y


def fwhm(y,x):
    y = y/amax(y)
    x_max = argmax(y)
    hwhm_r = -1
    hwhm_l = -1
    k = 0
    while hwhm_r == -1:
        k+=1
        if x_max+k>= x.size:
            break
        if y[x_max+k]<=0.5:
            hwhm_r = x[x_max+k]

    if hwhm_r == -1:
        hwhm_r= x[-1]
    else:
        hwhm_r = interp(0.5, y[x_max+k-3:x_max+k+3].flatten(),x[x_max+k-3:x_max+k+3].flatten())

    k = 0
    while hwhm_l == -1:
        k-=1
        if x_max+k<= 0:
            break
        if y[x_max+k]<=0.5:
            hwhm_l = x[x_max+k]

    if hwhm_l == -1:
        hwhm_l = x[0]
    else:
        hwhm_l = interp(0.5, y[x_max + k - 3:x_max + k + 3].flatten(), x[x_max + k - 3:x_max + k + 3].flatten())

    return (hwhm_r-hwhm_l)


def get_fwhm(image, x, y, coords_x, coords_y):
    f = interp2d(x, y, image, kind='linear')
    if coords_x == ':':
        coords_x = x
        cut = f(coords_x, coords_y)
        return fwhm(cut,coords_x)
    if coords_y == ':':
        coords_y = y
        cut = f(coords_x, coords_y)
        return fwhm(cut,coords_y)


if __name__ == '__main__':

    #############################################################################
    #############################################################################
    # Logging all prints

    # Get time stamp
    start0 = time.time()
    dt = datetime.datetime.fromtimestamp(start0).strftime('%Y-%m-%d_%H:%M:%S')

    # Initializing logging
    log = logging.getLogger('')
    log.setLevel(logging.INFO)
    format = logging.Formatter('%(levelname)s: %(message)s')

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(format)
    log.addHandler(ch)

    fh = logging.handlers.RotatingFileHandler(
        os.path.join(os.getcwd()) + '/AuxPlots_CSD_' + dt + '.log', maxBytes=(1048576 * 5), backupCount=7)
    fh.setFormatter(format)
    log.addHandler(fh)

    # ------------------------------------------------------------------------------------------------------------------
    # Intensity

    file = 'intensity7.0keV_15k_ME_intensity.dat'

    logging.info('>>> %s' % (file.split('/')[-1]))
    MI, mesh = srwl_uti_read_intens_ascii(file)
    logging.info('Intensity length: %d; nx: %d; ny: %d;'%(len(MI),mesh.nx,mesh.ny))

    MI = reshape(MI,(mesh.ny, mesh.nx))

    x = linspace(mesh.xStart,mesh.xFin, mesh.nx)
    y = linspace(mesh.yStart,mesh.yFin, mesh.ny)

    logging.info('Horizontal FWHM: %.2f [um]' % float(get_fwhm(MI, x, y, coords_x=':', coords_y=0)*1e6))
    logging.info('Vertical FWHM: %.2f [um]' % float(get_fwhm(MI, x, y, coords_x=0, coords_y=':')*1e6))

    image = b4pt.Image2Plot(MI/amax(MI), x*1e3, y*1e3)
    image.legends = ['', '(mm)', '(mm)']
    image.LaTex = True
    image.AspectRatio = True
    image.ColorScheme = 2
    image.plt_limits = [-0.05, 1.05]
    image.ax_limits = [-1,1,-1,1]
    image.FontsSizeScale = 1.3
    image.sort_class()
    b4pt.plot_2D_cuts(image, 'intensity_15kME_7keV.pdf', Enable=False, Silent=False)

    slc_Ih, axI_x = get_slice(MI, x, y, ':', 0)
    slc_Iv, axI_y = get_slice(MI, x, y, 0, ':')
    # ------------------------------------------------------------------------------------------------------------------
    # CSD_x

    file = 'CSDx7.0keV_15k_ME_intensity.dat'

    logging.info('>>> %s' % (file.split('/')[-1]))
    MI, mesh = srwl_uti_read_intens_ascii(file)
    logging.info('CSD: %d; nx: %d; ny: %d;'%(len(MI),mesh.nx,mesh.ny))

    CSD_Re = reshape(MI[0::2], (mesh.nx, mesh.nx))
    CSD_Im = reshape(MI[1::2], (mesh.nx, mesh.nx))

    CSD = CSD_Re + 1j*CSD_Im

    x = linspace(mesh.xStart,mesh.xFin, mesh.nx)
    y = linspace(mesh.xStart,mesh.xFin, mesh.nx)

    image = b4pt.Image2Plot(absolute(CSD)/amax(absolute(CSD)), x*1e3, y*1e3)
    image.legends = ['', '$x_1$ (mm)', '$x_2$ (mm)']
    image.LaTex = True
    image.AspectRatio = True
    image.ColorScheme = 7
    image.plt_limits = [-0.05, 1.05]
    image.ax_limits = [-1,1,-1,1]
    image.FontsSizeScale = 1.3
    image.sort_class()
    b4pt.plot_2D_cuts(image, 'CSDx_15kME_7keV.pdf', Enable=False, Silent=False)

    StokesParam = SRWLStokes(_arS=MI, _typeStokes='f', _eStart=mesh.eStart, _eFin=mesh.eFin, _ne=mesh.ne,
                             _xStart=mesh.xStart, _xFin=mesh.xFin, _nx=mesh.nx,
                             _yStart=mesh.yStart, _yFin=mesh.yFin, _ny=mesh.ny, _mutual=1)

    DoC = StokesParam.to_deg_coh()
    logging.info('DoC: %d; nx: %d; ny: %d;'%(len(DoC),mesh.nx,mesh.nx))
    DoC = reshape(DoC,(mesh.nx,mesh.nx))

    # logging.info('Horizontal FWHM: %.2f [um]' % (get_fwhm(DoC, x, y, coords_x=':', coords_y=0)*1e6))
    logging.info('Coherence length: %.2f [um]' % (get_fwhm(DoC, x, y, coords_x=0, coords_y=':')*1e6))

    image = b4pt.Image2Plot(absolute(DoC)/amax(DoC), x*1e3, y*1e3)
    image.legends = ['', '$(x_1+x_2)/2$ (mm)', '$(x_1-x_2)/2$ (mm)']
    image.LaTex = True
    image.AspectRatio = True
    image.ColorScheme = 5
    image.plt_limits = [-0.05, 1.05]
    image.ax_limits = [-1,1,-1,1]
    image.FontsSizeScale = 1.3
    image.sort_class()
    b4pt.plot_2D_cuts(image, 'DoCx_15kME_7keV.pdf', Enable=False, Silent=False)

    slc_DoCh, axDoC_x = get_slice(DoC, x, y, 0, ':')

    # ------------------------------------------------------------------------------------------------------------------
    # CSD_y

    file = 'CSDy7.0keV_15k_ME_intensity.dat'

    logging.info('>>> %s' % (file.split('/')[-1]))
    MI, mesh = srwl_uti_read_intens_ascii(file)
    logging.info('CSD: %d; nx: %d; ny: %d;'%(len(MI),mesh.nx,mesh.ny))

    CSD_Re = reshape(MI[0::2], (mesh.ny, mesh.ny))
    CSD_Im = reshape(MI[1::2], (mesh.ny, mesh.ny))

    CSD = CSD_Re + 1j*CSD_Im

    x = linspace(mesh.yStart,mesh.yFin, mesh.ny)
    y = linspace(mesh.yStart,mesh.yFin, mesh.ny)

    image = b4pt.Image2Plot(absolute(CSD)/amax(absolute(CSD)), x*1e3, y*1e3)
    image.legends = ['', '$y_1$ (mm)', '$y_2$ (mm)']
    image.LaTex = True
    image.AspectRatio = True
    image.ColorScheme = 7
    image.plt_limits = [-0.05, 1.05]
    image.ax_limits = [-1,1,-1,1]
    image.FontsSizeScale = 1.3
    image.sort_class()
    b4pt.plot_2D_cuts(image, 'CSDy_15kME_7keV.pdf', Enable=False, Silent=False)

    StokesParam = SRWLStokes(_arS=MI, _typeStokes='f', _eStart=mesh.eStart, _eFin=mesh.eFin, _ne=mesh.ne,
                             _xStart=mesh.xStart, _xFin=mesh.xFin, _nx=mesh.nx,
                             _yStart=mesh.yStart, _yFin=mesh.yFin, _ny=mesh.ny, _mutual=1)

    DoC = StokesParam.to_deg_coh()
    logging.info('DoC: %d; nx: %d; ny: %d;'%(len(DoC),mesh.ny,mesh.ny))
    DoC = reshape(DoC,(mesh.ny,mesh.ny))

    # logging.info('Horizontal FWHM: %.2f [um]' % (get_fwhm(DoC, x, y, coords_x=':', coords_y=0)*1e6))
    logging.info('Coherence length: %.2f [um]' % (get_fwhm(DoC, x, y, coords_x=0, coords_y=':')*1e6))

    image = b4pt.Image2Plot(absolute(DoC)/amax(absolute(DoC)), x*1e3, y*1e3)
    image.legends = ['', '$(y_1+y_2)/2$ (mm)', '$(y_1-y_2)/2$ (mm)']
    image.LaTex = True
    image.AspectRatio = True
    image.ColorScheme = 5
    image.plt_limits = [-0.05, 1.05]
    image.ax_limits = [-1,1,-1,1]
    image.FontsSizeScale = 1.3
    image.sort_class()
    b4pt.plot_2D_cuts(image, 'DoCy_15kME_7keV.pdf', Enable=False, Silent=False)

    slc_DoCv, axDoC_y = get_slice(DoC, x, y, 0, ':')

    # ------------------------------------------------------------------------------------------------------------------
    # 1D comp

    image = b4pt.Image2Plot(slc_Ih/amax(slc_Ih), axI_x*1e3)
    image.LaTex = True
    image.grid = True
    image.legends = ['', 'x (mm)', '(a.u.)']
    image.FontsSizeScale = 1.3
    image.ColorScheme = 1
    image.LineStyle = '-'
    image.label = 'intensity'
    image.LabelPos = 2
    image.ax_limits = [-1, 1, -0.05, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    image = b4pt.Image2Plot(slc_DoCh, axDoC_x*1e3)
    image.LaTex = True
    image.grid = True
    image.legends = ['', 'x (mm)', '(a.u.)']
    image.FontsSizeScale = 1.3
    image.ColorScheme = 2
    image.LineStyle = '-'
    image.label = 'DoC$_x$'
    image.LabelPos = 1
    image.FillBetween = True
    image.FillBetweenValue = 0
    image.alpha = 0.2
    image.ax_limits = [-1, 1, -0.05, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, 'Ih_vsDoCx_15kME_7keV.pdf', Enable=False, Hold=True)

    image = b4pt.Image2Plot(slc_Iv/amax(slc_Iv), axI_y*1e3)
    image.LaTex = True
    image.grid = True
    image.legends = ['', 'x (mm)', '(a.u.)']
    image.FontsSizeScale = 1.3
    image.ColorScheme = 1
    image.LineStyle = '-'
    image.label = 'intensity'
    image.LabelPos = 2
    image.ax_limits = [-1, 1, -0.05, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    image = b4pt.Image2Plot(slc_DoCv, axDoC_y*1e3)
    image.LaTex = True
    image.grid = True
    image.legends = ['', 'y (mm)', '(a.u.)']
    image.FontsSizeScale = 1.3
    image.ColorScheme = 2
    image.LineStyle = '-'
    image.label = 'DoC$_y$'
    image.LabelPos = 1
    image.FillBetween = True
    image.FillBetweenValue = 0
    image.alpha = 0.2
    image.ax_limits = [-1, 1, -0.05, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, 'Iv_vsDoCy_15kME_7keV.pdf', Enable=True, Hold=True)

