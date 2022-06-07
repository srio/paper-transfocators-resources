import numpy
from srxraylib.plot.gol import plot, set_qt
import matplotlib.pylab as plt
from scipy.signal import savgol_filter

set_qt()

# subdirectory = "/scisoft/users/srio/COMSYL-SLURM/id18_onelens/zeroemittance" # "Data7keV_200um"
subdirectory = "./zeroemittance/" # "Data7keV_200um"

beam_dimension_at_slit_in_um = 565  # needed for calculating Fresnel number
beam_dimension_at_source = 9.13 # needed for magnification and theoretical sizes
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
size_at_aperture = 565e-6
# a = 234e-6 / 2
# p = 65.0
# print("N = ",a**2 / (wavelength * p))

DISTANCE   = []
FWHM       = []
FWHM_SMOOTH       = []
ICENTER    = []
SIGMASF    = []
LEGEND     = []
FWHM_THEORY_SOURCE       = []
FWHM_THEORY_SLIT         = []
FWHM_THEORY_AVERAGE         = []
FWHM_RATIO_SOURCE         = []
FWHM_RATIO_SLIT         = []
APERTURE = []
BEST_FOCUS_POSITION = numpy.zeros(len(FACTOR))
BEST_FOCUS_FWHM = numpy.zeros(len(FACTOR))
BEST_FOCUS_I = numpy.zeros(len(FACTOR), dtype=int)
N2 = []
INTENSITY = []

for i in range(len(FACTOR)):
        filein = "%s/aperture_factor_%g.dat" % (subdirectory, FACTOR[i])
        print(">>>>> ", filein)
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; intensity1 = a1[:,2] ; icenter1 = a1[:,3]
        DISTANCE.append(distance1)
        FWHM.append(fwhm1)
        FWHM_SMOOTH.append( savgol_filter(savgol_filter(fwhm1, 7, 1), 7, 1) )
        ICENTER.append(icenter1)
        SIGMASF.append(float(FACTOR[i]))
        INTENSITY.append(intensity1[0])

        slit_size_in_um = beam_dimension_at_slit_in_um * float(FACTOR[i])
        s = slit_size_in_um * 1e-6 / 2
        p = 65.0
        pp = 36
        pa = p - pp
        # N2 =  p * s ** 2 / (wavelength * pp ** 2)
        N2.append( p * s ** 2 / (wavelength * pa ** 2) )
        print("Effect for slit aperture less than [um] = ", 2e6 * numpy.sqrt(pp ** 2 * wavelength / p))

        # magnification, theoretical sizes
        fwhm_theory_source = beam_dimension_at_source * distance1 / p
        fwhm_theory_slit = slit_size_in_um * distance1 / pa
        FWHM_RATIO_SOURCE.append(numpy.abs((fwhm1 - fwhm_theory_source) / fwhm1))
        FWHM_RATIO_SLIT.append(numpy.abs((fwhm1 - fwhm_theory_slit) / fwhm1))
        w1 = FACTOR[i]
        if w1 > 1:
                w1 = 1
        w1 = (w1) ** (1/15)
        w2 = 1.0 - w1
        wtot = w1 + w2
        w1 /= wtot
        w2 /= wtot
        fwhm_theory_average = fwhm_theory_source * (w1) + fwhm_theory_slit * (w2)
        # fwhm_theory_average = fwhm_theory_source * (N2/35) + fwhm_theory_slit * (1.0 - (N2/35))

        FWHM_THEORY_SOURCE.append(fwhm_theory_source)
        FWHM_THEORY_SLIT.append(fwhm_theory_slit)
        FWHM_THEORY_AVERAGE.append(fwhm_theory_average)


        APERTURE.append(size_at_aperture * FACTOR[i])

        i_best_focus = numpy.argmin(FWHM_SMOOTH[i])
        BEST_FOCUS_I[i] = i_best_focus
        BEST_FOCUS_FWHM[i] = FWHM_SMOOTH[i][i_best_focus]
        BEST_FOCUS_POSITION[i] = distance1[i_best_focus]

for i in range(len(FACTOR)):
        # LEGEND.append(r'$a$=%5.3g $T$=%5.3g; $n$=%5.3g; $N_F$=%4.3f' % (1e6*APERTURE[i],
        #                                                                 numpy.round(INTENSITY[i] / INTENSITY[-1], 2),
        #                                                                 FACTOR[i],
        #                                                                 N2[i]))

        # LEGEND.append(r'$a$=%5.3g $T$=%5.3g; $N_F$=%4.3f' % (1e6*APERTURE[i],
        #                                                                 numpy.round(INTENSITY[i] / INTENSITY[-1], 2),
        #                                                                 N2[i]))

        LEGEND.append(r'$a$=%5.3g $\mu$m; $T$=%5.3g' % (1e6*APERTURE[i], numpy.round(INTENSITY[i] / INTENSITY[-1], 2),))

show_fwhm =  True

if show_fwhm:
    fig, ax = plot(
            DISTANCE[0], FWHM_SMOOTH[0], # savgol_filter(savgol_filter(FWHM[0], 7, 1), 7, 1),
            DISTANCE[1], FWHM_SMOOTH[1], # savgol_filter(savgol_filter(FWHM[1], 7, 1), 7, 1),
            DISTANCE[2], FWHM_SMOOTH[2], # savgol_filter(savgol_filter(FWHM[2], 7, 1), 7, 1),
            DISTANCE[3], FWHM_SMOOTH[3], # savgol_filter(savgol_filter(FWHM[3], 7, 1), 7, 1),
            DISTANCE[4], FWHM_SMOOTH[4], # savgol_filter(savgol_filter(FWHM[4], 7, 1), 7, 1),
            DISTANCE[5], FWHM_SMOOTH[5], # savgol_filter(savgol_filter(FWHM[5], 7, 1), 7, 1),
            DISTANCE[6], FWHM_SMOOTH[6], # savgol_filter(savgol_filter(FWHM[6], 7, 1), 7, 1),
            DISTANCE[7], FWHM_SMOOTH[7], # savgol_filter(savgol_filter(FWHM[7], 7, 1), 7, 1),
            xrange=[8,62],
            ytitle=r'FWHM [$\mu$m]',yrange=[1,400],
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
# ax.set_title("Be Lens radius: %d um" % (lens_radius * 1e6))

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
ax.plot(BEST_FOCUS_POSITION, BEST_FOCUS_FWHM, color='black', linestyle='dashed')

# file_png = "oneTFund_200um.eps" #"scan_%s.pdf" % mirror_name
file_png = "oneTF_UndSource_RectSlit_R200um.eps" #"scan_%s.pdf" % mirror_name
plt.savefig(file_png)
print("File written to disk: %s" % file_png)


plt.show()
print(">>>", BEST_FOCUS_I[0] )
print(">>>", DISTANCE[0][BEST_FOCUS_I[0]])
print(">>>", FWHM_THEORY_SLIT[0][BEST_FOCUS_I[0]] )
# xx = DISTANCE[0][BEST_FOCUS_I[0]]
# yy = FWHM_THEORY_SLIT[0][BEST_FOCUS_I[0]]
# plot(BEST_FOCUS_POSITION, BEST_FOCUS_FWHM,
#         DISTANCE[0][BEST_FOCUS_I[0]], FWHM_THEORY_SLIT[0][BEST_FOCUS_I[0]],
#         DISTANCE[1][BEST_FOCUS_I[1]], FWHM_THEORY_SLIT[1][BEST_FOCUS_I[1]],
#         DISTANCE[2][BEST_FOCUS_I[2]], FWHM_THEORY_SLIT[2][BEST_FOCUS_I[2]],
#         DISTANCE[3][BEST_FOCUS_I[3]], FWHM_THEORY_SLIT[3][BEST_FOCUS_I[3]],
#         DISTANCE[4][BEST_FOCUS_I[4]], FWHM_THEORY_SLIT[4][BEST_FOCUS_I[4]],
#         DISTANCE[5][BEST_FOCUS_I[5]], FWHM_THEORY_SLIT[5][BEST_FOCUS_I[5]],
#         DISTANCE[6][BEST_FOCUS_I[6]], FWHM_THEORY_SLIT[6][BEST_FOCUS_I[6]],
#         DISTANCE[0][BEST_FOCUS_I[0]], FWHM_THEORY_SOURCE[0][BEST_FOCUS_I[0]],
#         DISTANCE[1][BEST_FOCUS_I[1]], FWHM_THEORY_SOURCE[1][BEST_FOCUS_I[1]],
#         DISTANCE[2][BEST_FOCUS_I[2]], FWHM_THEORY_SOURCE[2][BEST_FOCUS_I[2]],
#         DISTANCE[3][BEST_FOCUS_I[3]], FWHM_THEORY_SOURCE[3][BEST_FOCUS_I[3]],
#         DISTANCE[4][BEST_FOCUS_I[4]], FWHM_THEORY_SOURCE[4][BEST_FOCUS_I[4]],
#         DISTANCE[5][BEST_FOCUS_I[5]], FWHM_THEORY_SOURCE[5][BEST_FOCUS_I[5]],
#         DISTANCE[6][BEST_FOCUS_I[6]], FWHM_THEORY_SOURCE[6][BEST_FOCUS_I[6]],
#         marker=[None, 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'x', 'x', 'x', 'x', 'x', 'x', 'x'],
#         xrange=[15,30], )

