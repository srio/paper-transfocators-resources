
from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from srxraylib.plot.gol import plot

from oasys.util.oasys_util import get_fwhm
#
# create input_wavefront
#
#

input_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012,x_max=0.00012,number_of_points=10000)
input_wavefront.set_photon_energy(10000)

direction = "V"
if direction == "H":
    nmax = 49
    sigmaI = 3.03783e-05
    beta = 0.0922395
elif direction == "V":
    nmax = 4
    sigmaI = 5.00125e-06
    beta = 1.15083

for i in range(nmax):
    # input_wavefront.set_gaussian_hermite_mode(3.03783e-05, i, amplitude=1.0, shift=0.0, beta=0.0922395)
    # input_wavefront.set_gaussian_hermite_mode(5.00125e-06, i, amplitude=1.0, shift=0.0, beta=1.15083)
    input_wavefront.set_gaussian_hermite_mode(sigmaI, i, amplitude=1.0, shift=0.0, beta=beta)
    intensity = input_wavefront.get_intensity()
    if i == 0:
        intensities = intensity
    else:
        intensities += intensity

fwhm, _, _ = get_fwhm(intensities, input_wavefront.get_abscissas())
plot(1e6 * input_wavefront.get_abscissas(),intensities, title="FWHM = %f " % (1e6 * fwhm))