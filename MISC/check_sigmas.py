import numpy
import scipy.constants as codata


def propagate(fwhm, distance=36.0, wavelength=1e-10):
    sigma = fwhm / 2.355
    sigmap = wavelength / (4 * numpy.pi) / sigma


    sigma2 = sigmap * distance
    fwhm2 = sigma2 * 2.355
    print("Source: %g Propagated: %g um" % (1e6 * fwhm, 1e6 * fwhm2))





photon_energy = 7000
wavelength = codata.h * codata.c / codata.e / photon_energy
print("Wavelenth % g A" % (wavelength * 1e10))



fwhmh1 = 17.72e-6 # 17.52e-6
fwhmh2 = propagate(fwhmh1, wavelength=wavelength)

fwhmv1 = 10.78e-6 # 9.39e-6
fwhmv2 = propagate(fwhmv1, wavelength=wavelength)