import numpy
from srxraylib.plot.gol import plot
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
# theoretical
#
if False:
    wavelength = 1.77e-10
    D = 10.0
    separation = 20e-6
    width = 0.8e-6

    npoints = 2000
    x = numpy.linspace(-0.0015, 0.0015, npoints)

    arg = numpy.pi * width * x / (wavelength * D)
    g = 0.85
    I0 = (numpy.sinc(arg))**2
    I =  I0 * (1+ g * numpy.cos(2*numpy.pi/wavelength/D*separation*x))


    Iav = get_Iav(I)


    plot(1e6*x,I0,
         1e6*x,I,
         1e6*x,Iav,
        legend=["I0", "I", "Iav"])

    gfit = get_fitted_DoC(I)

    print("Goriginal, Gfitted = ", g, gfit)


#
# experimental
#
b = numpy.loadtxt("/home/srio/Oasys/doubleslit_spectral_density.dat")
II = b[:,1]
print("II Gfitted = ", get_fitted_DoC(II, do_plot=0))

