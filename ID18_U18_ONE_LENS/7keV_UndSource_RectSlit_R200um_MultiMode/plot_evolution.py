import numpy
from srxraylib.plot.gol import plot, set_qt
import matplotlib.pylab as plt

set_qt()
subdirectory = "."
lens_radius = 200e-6 # for plot title
direction = 'h'
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

for direction in ['h','v']:
    if direction == 'h':
        APERTURE = [40.3e-6, 85.1e-6, 145e-6, 1000e-6, -40.3e-6, -85.1e-6, -145e-6, -1000e-6]
        CF       = [90,      70,      50,     12]
    else:
        # APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6, 5000e-6, -25e-6, -227.0e-6, -506.7e-6, -1500e-6, -5000e-6]
        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6, -25e-6, -227.0e-6, -506.7e-6, -1500e-6]
        CF       = [99.9,        90,       70,      58]

    FILES = []
    for aperture in APERTURE:
            # src1, wf = main(aperture=aperture, distance=18.4168, number_of_points=number_of_points)

            FILES.append("aperture_%s_%g.dat" % (direction, 1e6 * aperture))



    import scipy.constants as codata
    wavelength = codata.h*codata.c/codata.e/7000


    DISTANCE   = []
    FWHM       = []
    ICENTER    = []
    SIGMASF    = []
    LEGEND     = []


    for i in range(len(APERTURE)):
            filein = FILES[i]
            print(">>>>> ", filein)
            a1 = numpy.loadtxt(filein)
            print(a1.shape)
            distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
            DISTANCE.append(distance1)
            FWHM.append(fwhm1)
            ICENTER.append(icenter1)

            if APERTURE[i] < 0:
                LEGEND.append(r'$a$=%4.1f (1st mode)' % (-1E6* APERTURE[i]))
            else:
                LEGEND.append(r'$a$=%4.1f (CF %g%s)' % (1E6 * APERTURE[i], CF[i], "%"))

    show_fwhm =  True

    if show_fwhm:
        fig, ax = plot(
                DISTANCE[0], 1e6*FWHM[0],
                DISTANCE[1], 1e6*FWHM[1],
                DISTANCE[2], 1e6*FWHM[2],
                DISTANCE[3], 1e6*FWHM[3],
                DISTANCE[4], 1e6*FWHM[4],
                DISTANCE[5], 1e6*FWHM[5],
                DISTANCE[6], 1e6*FWHM[6],
                DISTANCE[7], 1e6*FWHM[7],
                xrange=[8,62],
                ytitle='FWHM [$\mu$m]',
                xtitle="Distance from lens [m]",
                figsize=(15, 4),
                show=0,
                legend=LEGEND,
                color=['r','b','g','k','r','b','g','k'],
                linestyle=[None,None,None,None,'--','--','--','--'],
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

    ax.xaxis.grid()
    ax.yaxis.grid()

    file_png = "oneTF_UndSource_RectSlit_R200um_PartialCoherence_%s.eps" % (direction)
    plt.savefig(file_png)
    print("File written to disk: %s" % file_png)


    plt.show()

