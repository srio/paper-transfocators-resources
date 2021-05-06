import numpy
from gol import plot
import matplotlib.pylab as plt
import h5py



subdirectory = "." # "/home/esrf/reyesher/OASYS/EBSL1-ID18/focusing_paper"
beam_dimension_at_slit_in_um = 565  # needed for calculating Fresnel number
lens_radius = 200e-6 # for plot title

#
#
#
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'figure.autolayout': True})

params = {'legend.fontsize': 15,
          'legend.frameon': False,
          'legend.handlelength': 2,
          # 'axes.titlesize' : 24,
          'axes.labelsize' :   24,
          'lines.linewidth' :  3,
          'lines.markersize' : 10,
          'xtick.labelsize' :  25,
          'ytick.labelsize' :  25,
          # 'grid.color':       'r',
          # 'grid.linestyle':   '-',
          # 'grid.linewidth':     2,
          }
plt.rcParams.update(params)





FACTOR = [0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5] #, 2, 4, 6]

import scipy.constants as codata
wavelength = codata.h*codata.c/codata.e/7000
# a = 234e-6 / 2
# p = 65.0
# print("N = ",a**2 / (wavelength * p))

DISTANCE   = []
FWHM       = []
ICENTER    = []
SIGMASF    = []
LEGEND     = []

for i in range(len(FACTOR)):
    filein = "%s/rt_aperture_factor_%g.h5" % (subdirectory, FACTOR[i])
    print(">>>>> ", filein)
    a1 = h5py.File(filein,'r')
    #print(a1.shape)
    distance1 = numpy.array(a1['/caustic/fwhm/x'])
    fwhm1 = numpy.array(a1['/caustic/fwhm/y'])
    icenter1 = numpy.array(a1['/caustic/I0/y'])
    DISTANCE.append(distance1)
    FWHM.append(fwhm1)
    ICENTER.append(icenter1)
    SIGMASF.append(float(FACTOR[i]))
    slit_size_in_um = beam_dimension_at_slit_in_um * float(FACTOR[i])
    s = slit_size_in_um * 1e-6 / 2
    p = 65.0
    pp = 36
    pa = p - pp
    # N2 =  p * s ** 2 / (wavelength * pp ** 2)
    N2 = p * s ** 2 / (wavelength * pa ** 2)
    print("Effect for slit aperture less than [um] = ", 2e6 * numpy.sqrt(pp ** 2 * wavelength / p))
    LEGEND.append(r'$n$=%5.3g; N=%4.3f' % (FACTOR[i], N2))
    # LEGEND.append(r'$a$=%s a$_{FWHM}$' % (SIGMAS[i]))

#for i in range(len(FACTOR)):
#
#    plt.plot(DISTANCE[i],FWHM[i], label = LEGEND[i])
#
#plt.xlabel('Distance from lens [m]')
#plt.ylabel('FWHM [$\mu$m]')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = numpy.flipud(numpy.array(range(0,8,1)))
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper right')
##plt.legend()
#plt.show()
show_fwhm =  True

if show_fwhm:
    fig, ax = plot(
            DISTANCE[0], FWHM[0],
            DISTANCE[1], FWHM[1],
            DISTANCE[2], FWHM[2],
            DISTANCE[3], FWHM[3],
            DISTANCE[4], FWHM[4],
            DISTANCE[5], FWHM[5],
            DISTANCE[6], FWHM[6],
            DISTANCE[7], FWHM[7],
            ytitle="FWHM [$\mu$m]",
            xtitle="Distance from lens [m]",
            figsize=(15, 4),
            show=0,
            legend=LEGEND,
            ylog=1,)

else:
    fig, ax = plot(   DISTANCE[0], ICENTER[0],
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
            show=0,
            legend=LEGEND,
            ylog=1)

# plt.grid()
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_title("Be Lens radius: %d um" % (lens_radius * 1e6))

# F (source): 14.3507, p1: 65, q1: 18.4168
# F (slit): 14.3507, p1: 29, q1: 28.4089

# vertical lines
energy_keV = 7
R = lens_radius
import xraylib
xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", energy_keV, 1.85)).real
F = R / (2 * xrl_delta)
qsrc = 1/(1/F-1/p)
print("F (source): %g, p1: %g, q1: %g" % (F,p,qsrc))
p=65 -36.0
qslt = 1/(1/F-1/p)
print("F (slit): %g, p1: %g, q1: %g" % (F,p,qslt))
print("R_Be [mm]= ", 1e3*R)

ax.plot([qsrc,qsrc], [1e-9,10000])
ax.plot([qslt,qslt], [1e-9,10000])

file_png = "ray-tracing_FWHM_dist_from_lens.pdf"
plt.savefig(file_png)
print("File written to disk: %s" % file_png)

plt.show()
