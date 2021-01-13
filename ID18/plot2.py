import numpy
from srxraylib.plot.gol import plot

use_gaussian_slit = False
use_real_lens = True

SIGMAS = ["0.1","0.2","0.5","1.0","1.5","2.0","4.0","6.0"]

DISTANCE = []
FWHM = []
ICENTER = []
LEGEND = []

import scipy.constants as codata
wavelength = codata.h*codata.c/codata.e/10000
a = 234e-6 / 2
p = 65.0
print("N = ",a**2 / (wavelength * p))

if use_gaussian_slit:
        gauss_add = "G"
else:
        gauss_add = ""

if use_real_lens:
        real_lens_add = "R"
else:
        real_lens_add = ""

for i in range(len(SIGMAS)):
        filein = "tmp%s%s%s.dat" % (gauss_add, real_lens_add, SIGMAS[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        DISTANCE.append(distance1)
        FWHM.append(fwhm1)
        ICENTER.append(icenter1)

        slit_size_in_um = 125.0 / 2.35 * float(SIGMAS[i])
        s = slit_size_in_um * 1e-6 / 2
        pp = 35


        N2 =  p * s ** 2 / (wavelength * pp ** 2)
        print("Effect for slit aperture less than [um] = ", 2e6 * numpy.sqrt(pp ** 2 * wavelength / p))

        LEGEND.append(r'$a$=%s $\sigma_a$; N=%3.1f' % (SIGMAS[i], N2))

plot(   DISTANCE[0], FWHM[0],
        DISTANCE[1], FWHM[1],
        DISTANCE[2], FWHM[2],
        DISTANCE[3], FWHM[3],
        DISTANCE[4], FWHM[4],
        DISTANCE[5], FWHM[5],
        DISTANCE[6], FWHM[6],
        DISTANCE[7], FWHM[7],
        ytitle="FWHM [um]",
        xtitle="Distance from source [m]",
        figsize=(15, 4),
        show=0,
        legend=LEGEND,
        ylog=1)

print(len(DISTANCE),len(ICENTER),len(SIGMAS))
plot(   DISTANCE[0], ICENTER[0],
        DISTANCE[1], ICENTER[1],
        DISTANCE[2], ICENTER[2],
        DISTANCE[3], ICENTER[3],
        DISTANCE[4], ICENTER[4],
        DISTANCE[5], ICENTER[5],
        DISTANCE[6], ICENTER[6],
        DISTANCE[7], ICENTER[7],
        ytitle="Intensity on axis",
        xtitle="Distance from source [m]",
        figsize=(15, 4),
        show=1,
        legend=LEGEND,
        ylog=1)

# iMin1 = numpy.argmin(fwhm11)
# print("Minima found for: %g" % (distance1[iMin1]))
