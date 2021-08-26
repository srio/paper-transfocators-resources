#!/bin/python


import sys
sys.path.insert(0, './srw_python')
from srwlib import *

import numpy as np
import barc4plots.barc4plots as b4pt
from scipy.interpolate import interp2d


def read_srw_intensity_dat(file_name, bandwidth=1e-3, transmission=1, norm=False):

    image, mesh = srwl_uti_read_intens_ascii(file_name)
    image = np.reshape(image, (mesh.ny, mesh.nx))
    dx = (mesh.xFin - mesh.xStart)/mesh.nx * 1E3
    dy = (mesh.yFin - mesh.yStart)/mesh.ny * 1E3

    image = image*dx*dy*bandwidth*transmission/(1e-3)
    n = 1
    if norm is not False:
        if norm is True:
            image /= np.amax(image)
            n = np.amax(image)
        else:
            image /= norm
            n = norm

    x = np.linspace(mesh.xStart, mesh.xFin, mesh.nx)
    y = np.linspace(mesh.yStart, mesh.yFin, mesh.ny)

    return image, x, y, mesh, n


def read_srw_phase_dat(file_name, unwrap=False):

    image, mesh = srwl_uti_read_intens_ascii(file_name)
    image = np.reshape(image, (mesh.ny, mesh.nx))

    if unwrap:
        image = unwrap_phase(image, wrap_around=False)

    x = np.linspace(mesh.xStart, mesh.xFin, mesh.nx)
    y = np.linspace(mesh.yStart, mesh.yFin, mesh.ny)

    return image, x, y, mesh


def get_slice(_image, _x, _y, _coords_x, _coords_y):

    f = interp2d(_x, _y, _image, kind='linear')

    if _coords_x == ':':
        # coords_x = x
        _integrated = np.zeros(_x.size)
        _cut = f(_x, _coords_y)
        _axis = _x
        # return _cut, _x
    if _coords_y == ':':
        # _coords_y = y
        _integrated = np.zeros(_y.size)
        _cut = f(_coords_x, _y)
        _axis = _y

    for i in range(_integrated.size):
        _integrated[i] = _cut[i]
    return _integrated, _axis


def plts_analysis_c1():

    samples = np.asarray([1, 10, 50, 100, 500, 1000, 2000, 5000, 10000, 25000, 100000])
    pk_int = np.zeros(len(samples))
    err_rms = np.zeros(len(samples))
    err_rms_x = np.zeros(len(samples))
    err_rms_y = np.zeros(len(samples))

    lista = [
        './me_tests/id18_c01_7.0keV_0p001k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_0p010k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_0p050k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_0p100k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_0p500k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_1k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_2k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_5k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_10k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_25k_ME_intensity.dat',
        './me_tests/id18_c01_7.0keV_100k_ME_intensity.dat',
    ]

    img_ref, x, y, mesh, Imax = read_srw_intensity_dat('./me_tests/id18_c01_7.0keV_100k_ME_intensity.dat', norm=False)
    # I = np.amax(img_ref)
    cuty_ref, y = get_slice(img_ref, x, y, 0, ':')
    cutx_ref, x = get_slice(img_ref, x, y, ':', 0)
    k = 0
    for name in lista:
        print(name)
        img, x, y, mesh, Imax = read_srw_intensity_dat(name, norm=False)
        if k == 0:
            I = np.amax(img)
        pk_int[k] = np.amax(img/I)
        err_rms[k] = np.std(100*(img-img_ref)/img_ref)

        cuty, y = get_slice(img, x, y, 0, ':')
        cutx, x = get_slice(img, x, y, ':', 0)

        err_rms_x[k] = np.std(100*(cutx-cutx_ref)/cutx_ref)
        err_rms_y[k] = np.std(100*(cuty-cuty_ref)/cuty_ref)

        k+=1


    image = b4pt.Image2Plot(pk_int, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 0
    image.LineStyle = 'o-'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c01_7p0keV_peak_intensity.pdf', Enable=False, Hold=False)

    image = b4pt.Image2Plot(err_rms, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', 'std. dev. ($\%$)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 1
    image.LineStyle = 'o-'
    image.label = 'full image'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    image = b4pt.Image2Plot(err_rms_x, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', 'std. dev. ($\%$)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 2
    image.LineStyle = 'x-'
    image.label = 'ver. cut at $x = 0$'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    image = b4pt.Image2Plot(err_rms_y, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', 'std. dev. ($\%$)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 3
    image.label = 'hor. cut at $y = 0$'
    image.LineStyle = '^-'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c01_7p0keV_errs.pdf', Enable=True, Hold=True)


def plts_analysis_c3():


    samples = np.asarray([1, 10, 50, 100, 500, 1000, 2000, 5000, 10000, 25000, 100000])
    pk_int = np.zeros(len(samples))
    err_rms = np.zeros(len(samples))
    err_rms_x = np.zeros(len(samples))
    err_rms_y = np.zeros(len(samples))

    lista = [
        './me_tests/id18_c03_7.0keV_0p001k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_0p010k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_0p050k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_0p100k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_0p500k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_1k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_2k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_5k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_10k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_25k_ME_intensity.dat',
        './me_tests/id18_c03_7.0keV_100k_ME_intensity.dat',
    ]

    img_ref, x, y, mesh, Imax = read_srw_intensity_dat('./me_tests/id18_c03_7.0keV_100k_ME_intensity.dat', norm=False)
    # I = np.amax(img_ref)
    cuty_ref, y = get_slice(img_ref, x, y, 0, ':')
    cutx_ref, x = get_slice(img_ref, x, y, ':', 0)
    k = 0
    for name in lista:
        print(name)
        img, x, y, mesh, Imax = read_srw_intensity_dat(name, norm=False)
        if k == 0:
            I = np.amax(img)
        pk_int[k] = np.amax(img/I)
        err_rms[k] = np.std(100*(img-img_ref)/img_ref)

        cuty, y = get_slice(img, x, y, 0, ':')
        cutx, x = get_slice(img, x, y, ':', 0)

        err_rms_x[k] = np.std(100*(cutx-cutx_ref)/cutx_ref)
        err_rms_y[k] = np.std(100*(cuty-cuty_ref)/cuty_ref)

        k+=1


    image = b4pt.Image2Plot(pk_int, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 0
    image.LineStyle = 'o-'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c03_7p0keV_peak_intensity.pdf',Enable=False, Hold=False)

    image = b4pt.Image2Plot(err_rms, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', 'std. dev. ($\%$)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 1
    image.LineStyle = 'o-'
    image.label = 'full image'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    image = b4pt.Image2Plot(err_rms_x, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', 'std. dev. ($\%$)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 2
    image.LineStyle = 'x-'
    image.label = 'ver. cut at $x = 0$'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    image = b4pt.Image2Plot(err_rms_y, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', 'std. dev. ($\%$)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 3
    image.label = 'hor. cut at $y = 0$'
    image.LineStyle = '^-'
    image.Scale = 2
    image.ax_limits = [0.5, 200000, None, None]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c03_7p0keV_errs.pdf', Enable=True, Hold=True)

def plts_time():

    samples = np.asarray([100, 500, 1000, 2000, 5000, 10000, 25000, 100000])

    t_c1 = np.asarray([1.33333, 5, 9, 15, 39, 1*60+20, 3*60+21, 13*60+15])/60
    t_c2 = np.asarray([2.5, 9, 17, 30, 1*60+15, 2*60+35, 6*60+30, 4*(6*60+30)])/60

    image = b4pt.Image2Plot(t_c1-100, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 0
    image.LineStyle = 'o-'
    image.label = 'case 1'
    image.Scale = 2
    image.ax_limits = [50, 200000, -1, 28]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)
    image = b4pt.Image2Plot(t_c1, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 0
    image.LineStyle = 'o-'
    image.Scale = 2
    image.ax_limits = [50, 200000, -1, 28]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    image = b4pt.Image2Plot(t_c2, samples)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '(macro-electrons)', '(h)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 2
    image.LineStyle = '^-'
    image.label = 'case 3'
    image.Scale = 2
    image.ax_limits = [50, 200000, -1, 28]
    image.sort_class()
    b4pt.plot_1D(image, 'srw_time.pdf', Enable=True, Hold=True)



def plts_1dy():
    file_name = './me_tests/id18_c01_7.0keV_0p001k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    I = np.amax(cstx)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 101
    image.LineStyle = '-'
    image.label = '1'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    file_name = './me_tests/id18_c01_7.0keV_0p010k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 103
    image.LineStyle = '-'
    image.label = '10'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_0p050k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 105
    image.LineStyle = '-'
    image.label = '50'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_0p100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 107
    image.LineStyle = '-'
    image.label = '100'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_0p500k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 109
    image.LineStyle = '-'
    image.label = '500'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_1k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 111
    image.LineStyle = '-'
    image.label = '1k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_2k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 113
    image.LineStyle = '-'
    image.label = '2k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_5k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 115
    image.LineStyle = '-'
    image.label = '5k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_10k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 117
    image.LineStyle = '-'
    image.label = '10k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_25k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 119
    image.LineStyle = '-'
    image.label = '25k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = -1
    image.LineStyle = '-'
    image.label = '100k'
    image.LabelPos = 1
    image.ax_limits = [-25, 25, -0.02, 0.8]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c01_7p0keV_y_cut.pdf', Enable=False, Hold=True)

    # ------------------------------------------------------------------------------------------------------------------

    file_name = './me_tests/id18_c03_7.0keV_0p001k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    I = np.amax(cstx)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 101
    image.LineStyle = '-'
    image.label = '1'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    file_name = './me_tests/id18_c03_7.0keV_0p010k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 103
    image.LineStyle = '-'
    image.label = '10'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_0p050k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 105
    image.LineStyle = '-'
    image.label = '50'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_0p100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 107
    image.LineStyle = '-'
    image.label = '100'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_0p500k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 109
    image.LineStyle = '-'
    image.label = '500'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_1k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 111
    image.LineStyle = '-'
    image.label = '1k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_2k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 113
    image.LineStyle = '-'
    image.label = '2k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_5k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 115
    image.LineStyle = '-'
    image.label = '5k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_10k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 117
    image.LineStyle = '-'
    image.label = '10k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_25k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 119
    image.LineStyle = '-'
    image.label = '25k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, 0, ':')
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = -1
    image.LineStyle = '-'
    image.label = '100k'
    image.LabelPos = 1
    image.ax_limits = [-25, 25, -0.02, 0.5]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c03_7p0keV_y_cut.pdf', Enable=True, Hold=True)

def plts_1dx():
    file_name = './me_tests/id18_c01_7.0keV_0p001k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    I = np.amax(cstx)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 101
    image.LineStyle = '-'
    image.label = '1'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    file_name = './me_tests/id18_c01_7.0keV_0p010k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 103
    image.LineStyle = '-'
    image.label = '10'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_0p050k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 105
    image.LineStyle = '-'
    image.label = '50'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_0p100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 107
    image.LineStyle = '-'
    image.label = '100'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_0p500k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 109
    image.LineStyle = '-'
    image.label = '500'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_1k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 111
    image.LineStyle = '-'
    image.label = '1k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_2k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 113
    image.LineStyle = '-'
    image.label = '2k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_5k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 115
    image.LineStyle = '-'
    image.label = '5k'
    image.LabelPos = 1
    image.ax_limits = [-40, 40, None, None]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_10k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 117
    image.LineStyle = '-'
    image.label = '10k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_25k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 119
    image.LineStyle = '-'
    image.label = '25k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 0.8]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c01_7.0keV_100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = -1
    image.LineStyle = '-'
    image.label = '100k'
    image.LabelPos = 1
    image.ax_limits = [-35, 35, -0.02, 0.8]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c01_7.0keV_x_cut.pdf', Enable=False, Hold=True)

    # ------------------------------------------------------------------------------------------------------------------

    file_name = './me_tests/id18_c03_7.0keV_0p001k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    I = np.amax(cstx)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 101
    image.LineStyle = '-'
    image.label = '1'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=False)

    file_name = './me_tests/id18_c03_7.0keV_0p010k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 103
    image.LineStyle = '-'
    image.label = '10'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_0p050k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 105
    image.LineStyle = '-'
    image.label = '50'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_0p100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 107
    image.LineStyle = '-'
    image.label = '100'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_0p500k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 109
    image.LineStyle = '-'
    image.label = '500'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_1k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 111
    image.LineStyle = '-'
    image.label = '1k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_2k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 113
    image.LineStyle = '-'
    image.label = '2k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_5k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 115
    image.LineStyle = '-'
    image.label = '5k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_10k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 117
    image.LineStyle = '-'
    image.label = '10k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 1.05]
    image.sort_class()
    b4pt.plot_1D(image, Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_25k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = 119
    image.LineStyle = '-'
    image.label = '25k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 0.5]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c03_7.0keV_x_cut.pdf', Enable=False, Hold=True)

    file_name = './me_tests/id18_c03_7.0keV_100k_ME_intensity.dat'
    print(file_name)
    cstx, x, y, mesh, Imax = read_srw_intensity_dat(file_name, norm=False)
    cuty, axisy = get_slice(cstx/I, x, y, ':',0)
    image = b4pt.Image2Plot(cuty, axisy*1e6)
    image.LaTex = True
    image.grid = True
    image.legends = ['', '($\mu$m)', '(a.u.)']
    image.FontsSizeScale = 1.1
    image.ColorScheme = -1
    image.LineStyle = '-'
    image.label = '100k'
    image.LabelPos = 1
    image.ax_limits = [-75, 75, -0.02, 0.5]
    image.sort_class()
    b4pt.plot_1D(image, 'id18_c03_7.0keV_x_cut.pdf', Enable=True, Hold=True)

if __name__ == '__main__':
    #
    # plts_analysis_c1()
    # plts_analysis_c3()
    plts_time()
    # plts_1dy()
    # plts_1dx()

