import numpy
from srxraylib.plot.gol import plot, set_qt
import matplotlib.pylab as plt
from scipy.signal import savgol_filter
import scipy.constants as codata
import xraylib


def get_zero_emittance():
        # subdirectory = "/scisoft/users/srio/COMSYL-SLURM/id18_onelens/zeroemittance" # "Data7keV_200um"
        subdirectory = "./7keV_UndSource_RectSlit_R200um_Letter/zeroemittance/"  # "Data7keV_200um"

        beam_dimension_at_slit_in_um = 565  # needed for calculating Fresnel number
        beam_dimension_at_source = 9.13  # needed for magnification and theoretical sizes
        lens_radius = 200e-6  # for plot title

        FACTOR = [0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5]  # , 2, 4, 6]

        import scipy.constants as codata
        wavelength = codata.h * codata.c / codata.e / 7000
        size_at_aperture = 565e-6
        # a = 234e-6 / 2
        # p = 65.0
        # print("N = ",a**2 / (wavelength * p))

        DISTANCE = []
        FWHM = []
        FWHM_SMOOTH = []
        ICENTER = []
        SIGMASF = []
        LEGEND = []
        FWHM_THEORY_SOURCE = []
        FWHM_THEORY_SLIT = []
        FWHM_THEORY_AVERAGE = []
        FWHM_RATIO_SOURCE = []
        FWHM_RATIO_SLIT = []
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
                distance1 = a1[:, 0];
                fwhm1 = a1[:, 1];
                intensity1 = a1[:, 2];
                icenter1 = a1[:, 3]
                DISTANCE.append(distance1)
                FWHM.append(fwhm1)
                FWHM_SMOOTH.append(savgol_filter(savgol_filter(fwhm1, 7, 1), 7, 1))
                ICENTER.append(icenter1)
                SIGMASF.append(float(FACTOR[i]))
                INTENSITY.append(intensity1[0])

                slit_size_in_um = beam_dimension_at_slit_in_um * float(FACTOR[i])
                s = slit_size_in_um * 1e-6 / 2
                p = 65.0
                pp = 36
                pa = p - pp
                # N2 =  p * s ** 2 / (wavelength * pp ** 2)
                N2.append(p * s ** 2 / (wavelength * pa ** 2))
                print("Effect for slit aperture less than [um] = ", 2e6 * numpy.sqrt(pp ** 2 * wavelength / p))

                # magnification, theoretical sizes
                fwhm_theory_source = beam_dimension_at_source * distance1 / p
                fwhm_theory_slit = slit_size_in_um * distance1 / pa
                FWHM_RATIO_SOURCE.append(numpy.abs((fwhm1 - fwhm_theory_source) / fwhm1))
                FWHM_RATIO_SLIT.append(numpy.abs((fwhm1 - fwhm_theory_slit) / fwhm1))
                w1 = FACTOR[i]
                if w1 > 1:
                        w1 = 1
                w1 = (w1) ** (1 / 15)
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

                # LEGEND.append(r'$a$=%5.3g $\mu$m; $T$=%5.3g' % (1e6 * APERTURE[i], numpy.round(INTENSITY[i] / INTENSITY[-1], 2),))
                LEGEND.append(r'$a$=%5.3g ($T$=%5.3g)' % (1e6 * APERTURE[i], numpy.round(INTENSITY[i] / INTENSITY[-1], 2),))

        return DISTANCE, FWHM_SMOOTH, LEGEND, BEST_FOCUS_POSITION, BEST_FOCUS_FWHM


def get_emittances(direction='h'):
        if direction == 'h':
                APERTURE = [40.3e-6, 85.1e-6, 145e-6, 1000e-6, -40.3e-6, -85.1e-6, -145e-6, -1000e-6]
                CF = [90, 70, 50, 12]
        else:
                # APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6, 5000e-6, -25e-6, -227.0e-6, -506.7e-6, -1500e-6, -5000e-6]
                APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6, -25e-6, -227.0e-6, -506.7e-6, -1500e-6]
                CF = [99.9, 90, 70, 58]

        FILES = []
        for aperture in APERTURE:
                # src1, wf = main(aperture=aperture, distance=18.4168, number_of_points=number_of_points)

                FILES.append("./7keV_UndSource_RectSlit_R200um_MultiMode/aperture_%s_%g.dat" % (direction, 1e6 * aperture))

        import scipy.constants as codata
        wavelength = codata.h * codata.c / codata.e / 7000

        DISTANCE = []
        FWHM = []
        ICENTER = []
        SIGMASF = []
        LEGEND = []

        for i in range(len(APERTURE)):
                filein = FILES[i]
                print(">>>>> ", filein)
                a1 = numpy.loadtxt(filein)
                print(a1.shape)
                distance1 = a1[:, 0];
                fwhm1 = a1[:, 1];
                icenter1 = a1[:, 3]
                DISTANCE.append(distance1)
                FWHM.append(fwhm1)
                ICENTER.append(icenter1)

                # if APERTURE[i] < 0:
                #         LEGEND.append(r'$a$=%4.1f $\mu$m (1st mode)' % (-1E6 * APERTURE[i]))
                # else:
                #         LEGEND.append(r'$a$=%4.1f $\mu$m (CF %g%s)' % (1E6 * APERTURE[i], CF[i], "%"))

                if APERTURE[i] < 0:
                        LEGEND.append(r'$a$=%4.1f (1st mode)' % (-1E6 * APERTURE[i]))
                else:
                        LEGEND.append(r'$a$=%4.1f (CF %g%s)' % (1E6 * APERTURE[i], CF[i], "%"))

        return DISTANCE, FWHM, LEGEND

if __name__ == "__main__":

        title_fontsize = 34
        legend_fontsize = 18
        legend_position = (0.78, 1.)


        plt.rcParams.update(plt.rcParamsDefault)
        plt.rcParams.update({'figure.autolayout': True})
        plt.rcParams.update({'axes.linewidth' : 5})

        params = {'legend.fontsize': 18,
                  'legend.frameon': False,
                  'legend.handlelength': 2,
                  # 'axes.titlesize' : 24,
                  'axes.labelsize': 35,
                  'lines.linewidth': 3,
                  'lines.markersize': 10,
                  'xtick.labelsize': 35,
                  'ytick.labelsize': 35,
                  # 'grid.color':       'r',
                  # 'grid.linestyle':   '-',
                  # 'grid.linewidth':     2,
                  }
        plt.rcParams.update(params)

        DISTANCE, FWHM_SMOOTH, LEGEND, BEST_FOCUS_POSITION, BEST_FOCUS_FWHM = get_zero_emittance()

        fig = plt.figure(figsize=(23, 11))
        gs = fig.add_gridspec(3, hspace=0)
        axs = gs.subplots(sharex=True, sharey=True)

        # for i in range(3):
        #         axs[i].spines['right'].set_visible(False)
                # axs[i].spines['top'].set_visible(False)

        ##########################################################################################################
        # zero emittance
        ##########################################################################################################

        # fig, ax = plot(
        #         DISTANCE[0], FWHM_SMOOTH[0],  # savgol_filter(savgol_filter(FWHM[0], 7, 1), 7, 1),
        #         DISTANCE[1], FWHM_SMOOTH[1],  # savgol_filter(savgol_filter(FWHM[1], 7, 1), 7, 1),
        #         DISTANCE[2], FWHM_SMOOTH[2],  # savgol_filter(savgol_filter(FWHM[2], 7, 1), 7, 1),
        #         DISTANCE[3], FWHM_SMOOTH[3],  # savgol_filter(savgol_filter(FWHM[3], 7, 1), 7, 1),
        #         DISTANCE[4], FWHM_SMOOTH[4],  # savgol_filter(savgol_filter(FWHM[4], 7, 1), 7, 1),
        #         DISTANCE[5], FWHM_SMOOTH[5],  # savgol_filter(savgol_filter(FWHM[5], 7, 1), 7, 1),
        #         DISTANCE[6], FWHM_SMOOTH[6],  # savgol_filter(savgol_filter(FWHM[6], 7, 1), 7, 1),
        #         DISTANCE[7], FWHM_SMOOTH[7],  # savgol_filter(savgol_filter(FWHM[7], 7, 1), 7, 1),
        #         xrange=[8, 62],
        #         ytitle=r'FWHM [$\mu$m]', yrange=[1, 400],
        #         xtitle="Distance from lens [m]",
        #         figsize=(15, 4),
        #         show=0,
        #         legend=LEGEND,
        #         ylog=1, )

        for i in range(len(DISTANCE)):
                axs[0].plot(DISTANCE[i], FWHM_SMOOTH[i], label=LEGEND[i])

        # axs[0].plot(DISTANCE[1], FWHM_SMOOTH[1], )
        # axs[0].plot(DISTANCE[2], FWHM_SMOOTH[2], )
        # axs[0].plot(DISTANCE[3], FWHM_SMOOTH[3], )
        # axs[0].plot(DISTANCE[4], FWHM_SMOOTH[4], )
        # axs[0].plot(DISTANCE[5], FWHM_SMOOTH[5], )
        # axs[0].plot(DISTANCE[6], FWHM_SMOOTH[6], )
        # axs[0].plot(DISTANCE[7], FWHM_SMOOTH[7], )

        axs[0].set_yscale("log")
        plt.xlim([8, 62])
        plt.ylim([1, 400])


        #
        # vertical lines
        #
        energy_keV = 7
        lens_radius = 200e-6  # for plot title
        p = 65.0

        R = lens_radius

        xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", energy_keV, 1.85)).real
        F = R / (2 * xrl_delta)
        qsrc = 1/(1/F-1/p)
        print("F (source): %g, p1: %g, q1: %g" % (F,p,qsrc))
        p=65 -36.0
        qslt = 1/(1/F-1/p)
        print("F (slit): %g, p1: %g, q1: %g" % (F,p,qslt))
        print("R_Be [mm]= ", 1e3*R)


        axs[0].plot(BEST_FOCUS_POSITION, BEST_FOCUS_FWHM, color='black', linestyle='dashed')
        for i in range(3):
                axs[i].plot([qsrc,qsrc], [1e-9,10000], color='blue')
                axs[i].plot([qslt,qslt], [1e-9,10000], color='orange')



        axs[0].legend(bbox_to_anchor=legend_position, fontsize=legend_fontsize)

        plt.text(10, 300000, r"$a) \epsilon$=0", fontsize=36)

        ##########################################################################################################
        # H emittance
        ##########################################################################################################

        DISTANCE, FWHM, LEGEND = get_emittances(direction='h')

        for i in range(len(DISTANCE)):
                print(FWHM[i])
                axs[1].plot(DISTANCE[i], 1e6 * FWHM[i], label=LEGEND[i])

        axs[1].legend(bbox_to_anchor=legend_position, fontsize=legend_fontsize)


        plt.text(10, 1000, r"b) $\epsilon$=130 pm.rad", fontsize=36)

        ##########################################################################################################
        # V emittance
        ##########################################################################################################

        DISTANCE, FWHM, LEGEND = get_emittances(direction='v')

        for i in range(len(DISTANCE)):
                print(FWHM[i])
                axs[2].plot(DISTANCE[i], 1e6 * FWHM[i], label=LEGEND[i])

        axs[2].legend(bbox_to_anchor=legend_position, fontsize=legend_fontsize)

        ytitle=r'                                                 FWHM [$\mu$m]'
        xtitle="Distance from lens [m]"

        plt.xlabel(xtitle, fontsize=title_fontsize)
        plt.ylabel(ytitle, fontsize=title_fontsize)

        #
        # final
        #

        for i in range(3):
                axs[i].xaxis.grid()
                axs[i].yaxis.grid()

        # plt.text(10, 13000, r"$\epsilon$=130 pm.rad", fontsize=36)
        plt.text(10, 2, r"c) $\epsilon$=10 pm.rad", fontsize=36)




        file_png = "evolution.eps" #"scan_%s.pdf" % mirror_name
        plt.savefig(file_png)
        print("File written to disk: %s" % file_png)

        plt.show()
