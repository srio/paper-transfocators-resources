import numpy

from srxraylib.plot.gol import plot, plot_image
import matplotlib.pylab as plt
import xraylib
from silx.io.specfile import SpecFile


print("h/v     slit    F1        F2        R1         R2       FWHM     idx")

A0 = []
A1 = []
A1N = []
LEGEND = []

for direction in ['h','v']:
    # direction = 'v'


    if direction == 'h':
        APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]
        Selected_f1 = [46, 25]

        # aperture = APERTURE[0]
        # Cases = [1,2]

        aperture = APERTURE[1]
        Cases = [3,4]
    else:
        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]
        Selected_f1 = [15, 42]

        # aperture = APERTURE[1]
        # Cases = [1,2]

        aperture = APERTURE[2]
        Cases = [3,4]

    fileroot = "case7keV_%s" % direction
    subdirectory = "./sizes_slit%g_%s" % (1e6 * aperture, direction)


    a = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_slit%g_%s.dat" % (1e6 * aperture, direction))
    F1 = a[:, 0].copy()
    F2 = a[:, 1].copy()

    R1 = []
    R2 = []
    for i in range(F1.size):
        xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
        R1.append(F1[i] * (2 * xrl_delta))
        R2.append(F2[i] * (2 * xrl_delta))


    # plot(F1, F2, title="%s  slit: %g um" % (direction, aperture*1e6))




    for selected_f1 in Selected_f1:
        index = numpy.argmin(  numpy.abs(F1-selected_f1) )

        # try:
        if True:
            filename = "sizes_slit%g_%s/case7keV_%s_spectral_density_%03d.dat" % (
                aperture * 1e6, direction, direction, index)
            sf = SpecFile(filename)
            s1 = sf[0]
            fwhm = float(s1.scan_header_dict["UFWHM"])
            A0.append(s1.data[0,:])
            A1.append(s1.data[1,:])
            A1N.append(s1.data[1, :] / s1.data[1, :].max())
            LEGEND.append("%s a=%g um" % (direction, 1e6 * aperture))
        # except:
        #     fwhm = 0


        print("%s &      %4.1f & %3.1f &     %3.1f &     %3.1f &     %3.1f &     %3.1f &     %d \\\\" % (direction, 1e6*aperture, F1[index], F2[index],
                                                                                  1e6*R1[index], 1e6*R2[index], fwhm, index))



plot(A0[0], A1N[0],
     A0[2], A1N[2],
     legend=["case %d h" % Cases[0], "case %d v" % Cases[0]], figsize=[10,10], show=0)

filename = "case_%d_profiles.eps" % Cases[0]
plt.savefig(filename)
print("File written to disk: %s" % filename)
plt.show()


plot(A0[1], A1N[1],
     A0[3], A1N[3],
     legend=["case %d h" % Cases[1], "case %d v" % Cases[1]], figsize=[10,10], show=0)

filename = "case_%d_profiles.eps" % Cases[1]
plt.savefig(filename)
print("File written to disk: %s" % filename)
plt.show()


plot_image( numpy.outer(A1N[0], A1N[2]), A0[0], A0[2], title="case %d" % Cases[0], figsize=[10,10], aspect=None, show=0)
filename = "case_%d_image.eps" % Cases[0]
plt.savefig(filename)
print("File written to disk: %s" % filename)
plt.show()

plot_image( numpy.outer(A1N[1], A1N[3]), A0[1], A0[3], title="case %d" % Cases[1], figsize=[10,10], aspect=None, show=0)
filename = "case_%d_image.eps" % Cases[1]
plt.savefig(filename)
print("File written to disk: %s" % filename)
plt.show()

