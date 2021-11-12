import numpy

from srxraylib.plot.gol import plot, plot_image
import matplotlib.pylab as plt
import xraylib
from silx.io.specfile import SpecFile

import matplotlib
matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)
import pylab as plt
params = {'legend.fontsize': 20,
          'legend.handlelength': 2}
plt.rcParams.update(params)



print("h/v     slit    F1        F2        R1         R2 ")



Cases = ['1h','1v','2h','2v','3h','3v','4h','4v']
APERTURE = [40.3e-6,227.0e-6,40.3e-6,227.0e-6,85.1e-6, 506.7e-6,85.1e-6, 506.7e-6,]
Selected_f1 = [46.1, 15.0, 25.1, 42.2, 46.1, 85.2, 25.1, 42.2]
directory = "trajectories_precalculated"
directory = "trajectories_precalculated_refined"

for i,case in enumerate(Cases):

    # direction = 'v'


    direction = Cases[i][1]
    selected_f1 = Selected_f1[i]
    aperture = APERTURE[i]

    filename = "%s/f1_vs_f2_slit%g_%s.dat" % (directory, 1e6 * aperture, direction)
    a = numpy.loadtxt(filename)
    F1 = a[:, 0].copy()
    F2 = a[:, 1].copy()

    R1 = []
    R2 = []
    for i in range(F1.size):
        xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
        R1.append(F1[i] * (2 * xrl_delta))
        R2.append(F2[i] * (2 * xrl_delta))



    index = numpy.argmin(  numpy.abs(F1-selected_f1) )


    print("%s &      %4.1f & %3.1f &     %3.1f &     %3.1f &     %3.1f &     %d \\\\" % (case, 1e6*aperture, F1[index], F2[index],
                                                                              1e6*R1[index], 1e6*R2[index], index))



