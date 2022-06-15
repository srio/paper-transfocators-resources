import numpy
from srxraylib.plot.gol import plot, set_qt

set_qt()

import matplotlib.pylab as plt
from scipy.signal import savgol_filter


def get_f2(f1=28.2,
           position_source=0.0,
           position_lens1=65.0,
           position_lens2= 170.0,
           position_sample=200.0,
           verbose=True):

    p1 = position_lens1 - position_source
    q1 = 1 / (1 / f1 - 1 / p1)


    p2 = position_lens2 - (p1 + q1)
    q2 = position_sample - position_lens2

    f2 = 1.0 / (1 / p2 + 1 / q2)

    if verbose:
        D = position_lens2 - position_lens1
        print("D: %g, q1+p2: %g" % (D, q1+p2))
        print("p1: %g" % p1)
        print("q1: %g" % q1)
        print("p2: %g" % p2)
        print("q2: %g" % q2)
        print("D: %g, Sum: %g" % (D, p1+q1+p2+q2))

    M = (q1 / p1) * (q2 / p2)
    return f2, M

#
#
#


#
#
#
def plot_sizes_analytical_and_numerical(filename_root=None):


    sourcesize_h = 70.57e-6
    sourcesize_v = 15.02e-6




    F1 = numpy.linspace(5, 100, 500)

    F2theory1 = []
    F2theory2 = []
    Msource_at_id = []
    Msource_at_slit = []
    ff_source_at_id = []
    ff_source_at_slit = []

    for index, f1 in enumerate(F1):
        ff_id, mm_source_at_id = get_f2(f1=f1,
                                 position_source=0.0,  # source at source
                                 position_lens1=65.0,
                                 position_lens2=170.0,
                                 position_sample=200.0,
                                 verbose=False)

        ff_slit, mm_source_at_slit = get_f2(f1=F1[index],
                                 position_source=36.0, #0.0,  # source at slit
                                 position_lens1=65.0,
                                 position_lens2=170.0,
                                 position_sample=200.0,
                                 verbose=False)

        if ff_id < 0:
            ff_id = numpy.nan
            mm_source_at_id = numpy.nan

        if ff_slit < 0:
            ff_slit = numpy.nan
            mm_source_at_slit = numpy.nan

        Msource_at_id.append(mm_source_at_id)
        Msource_at_slit.append(mm_source_at_slit)

        F2theory1.append(ff_source_at_id)
        F2theory2.append(ff_source_at_slit)
        ff_source_at_id.append(ff_id)
        ff_source_at_slit.append(ff_slit)


    Msource_at_id   = numpy.abs(numpy.array(Msource_at_id))
    Msource_at_slit = numpy.abs(numpy.array(Msource_at_slit))


    if True:

        aperture_h = 40.3e-6 / ( numpy.sqrt(2*numpy.pi) / 2.355)  # corrected for rectangle with same area as Gaussian
        aperture_v = 25e-6   / ( numpy.sqrt(2*numpy.pi) / 2.355)  # corrected for rectangle with same area as Gaussian

        # aperture_h = 40.3e-6 * 0.7  # corrected a huevo
        # aperture_v = 25e-6   * 0.7  # corrected a huevo


        Hsrc = Msource_at_id * sourcesize_h * 1e6
        Hslt = Msource_at_slit * aperture_h * 1e6
        Vsrc = Msource_at_id * sourcesize_v * 1e6
        Vslt = Msource_at_slit * aperture_v * 1e6

        if False:
            Hsrc = savgol_filter(savgol_filter(Hsrc, 7, 1), 7, 1)
            Hslt = savgol_filter(savgol_filter(Hslt, 7, 1), 7, 1)
            Vsrc = savgol_filter(savgol_filter(Vsrc, 7, 1), 7, 1)
            Vslt = savgol_filter(savgol_filter(Vslt, 7, 1), 7, 1)



        APERTURE_H = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]

        APERTURE_V = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]

        HH0 = numpy.loadtxt("plot_combined_results_H_0.dat", skiprows=1)
        VV0 = numpy.loadtxt("plot_combined_results_V_0.dat", skiprows=1)
        # CF = [0.9, 0.999]
        # aperture_h = 40.3e-6
        # aperture_v = 25e-6


        HH1 = numpy.loadtxt("plot_combined_results_H_1.dat", skiprows=1)
        VV1 = numpy.loadtxt("plot_combined_results_V_1.dat", skiprows=1)
        # CF = [0.9, 0.9]
        # aperture_h = 40.3e-6
        # aperture_v = 227.0e-6
        #
        #
        HH2 = numpy.loadtxt("plot_combined_results_H_2.dat", skiprows=1)
        VV2 = numpy.loadtxt("plot_combined_results_V_2.dat", skiprows=1)
        CF = [0.7, 0.7]
        aperture_h = 85.1e-6
        aperture_v = 506.7e-6
        # #
        #
        HH3 = numpy.loadtxt("plot_combined_results_H_3.dat", skiprows=1)
        VV3 = numpy.loadtxt("plot_combined_results_V_3.dat", skiprows=1)
        # CF = [0.13, 0.58]
        # aperture_h = 1000e-6
        # aperture_v = 1500e-6


        # Hsrc = Msource_at_id * sourcesize_h * 1e6
        # Hslt = Msource_at_slit * aperture_h * 1e6
        # Vsrc = Msource_at_id * sourcesize_v * 1e6
        # Vslt = Msource_at_slit * aperture_v * 1e6
        #
        # # linear interpolation
        # m_H = (Hslt - Hsrc) / (1 - 0.13)
        # m_V = (Vslt - Vsrc) / (1 - 0.58)
        # y_H = m_H * (CF[0] - 0.13) + Hsrc
        # y_V = m_V * (CF[1] - 0.58) + Vsrc


        # h

        if True:
            fig2, ax2 = plot(
                HH0[:, 1][::5], HH0[:, 7][::5],
                HH1[:, 1][::5], HH1[:, 7][::5],
                HH2[:, 1][::5], HH2[:, 7][::5],
                HH3[:, 1][::5], HH3[:, 7][::5],
                numpy.array(F1), Hsrc,
                numpy.array(F1), Hslt,
                color=['b','cyan','pink','r','r','b'],
                 marker=['o','o','o','o',None,None],
                 linestyle=['','','','','--','--'],
                 # yrange=[1,1e3], ylog=1,
                 yrange=[0,50], ylog=0,
                 legend=["a=%g" % (1e6*APERTURE_H[0]),"a=%g" % (1e6*APERTURE_H[1]),"a=%g" % (1e6*APERTURE_H[2]),"a=%g" % (1e6*APERTURE_H[3]),
                         "Analytical source at Und", "Analytical source at Slit"],
                 xtitle="F1 [m]", ytitle="FWHM [$\mu$m]", title="",
                 figsize=(16, 8),
                 show=0)

            # ax2.fill_between(HH0[:, 1], HH0[:, -1], HH3[:, -1])
            ax2.xaxis.grid()
            ax2.yaxis.grid()

            if filename_root is not None:
                filename = "%s_h.eps" % (filename_root, )
                plt.savefig(filename)
                print("File written to disk: %s" % filename)

            plt.show()

        # v
        if True:
            fig2, ax2 = plot(
                VV0[:, 1][::3], VV0[:, 7][::3],
                VV1[:, 1][::3], VV1[:, 7][::3],
                VV2[:, 1][::3], VV2[:, 7][::3],
                VV3[:, 1][::3], VV3[:, 7][::3],
                numpy.array(F1)[::15], Vsrc[::15],
                numpy.array(F1)[::15], Vslt[::15],
                color=['b','cyan','pink','r','r','b'],
                 marker=['o','o','o','o',None,None],
                 linestyle=['','','','','--','--'],
                 # yrange=[1,1e3], ylog=1,
                 yrange=[0,50], ylog=0,
                 legend=["a=%g" % (1e6*APERTURE_V[0]),"a=%g" % (1e6*APERTURE_V[1]),"a=%g" % (1e6*APERTURE_V[2]),"a=%g" % (1e6*APERTURE_V[3]),
                         "Analytical source at Und", "Analytical source at Slit"],
                 xtitle="F1 [m]", ytitle="FWHM [$\mu$m]", title="",
                 figsize=(16,8),
                 show=0)

            # ax2.fill_between(VV0[:, 1], VV0[:, -1], VV3[:, -1])
            ax2.xaxis.grid()
            ax2.yaxis.grid()

            if filename_root is not None:
                filename = "%s_v.eps" % (filename_root)
                plt.savefig(filename)
                print("File written to disk: %s" % filename)

            plt.show()

            return VV0[:, 1][::3], VV0[:, 4][::3], VV0[:, 5][::3]
if __name__ == "__main__":

    plt.rcParams.update(plt.rcParamsDefault)
    # plt.rcParams.update({'figure.autolayout': True})
    # plt.rcParams.update({'axes.linewidth': 5})

    params = {'legend.fontsize': 20,
              'legend.frameon': False,
              'legend.handlelength': 2,
              # 'axes.titlesize' : 24,
              'axes.labelsize': 28,
              'lines.linewidth': 3,
              'lines.markersize': 10,
              'xtick.labelsize': 28,
              'ytick.labelsize': 28,
              # 'grid.color':       'r',
              # 'grid.linestyle':   '-',
              # 'grid.linewidth':     2,
              }
    plt.rcParams.update(params)

    f1, f2, fwhm = plot_sizes_analytical_and_numerical(filename_root="sizes")

    f = open('new_sizes0.dat','w')
    for i in range(f1.size):
        f.write("%g  %g  %g\n" % (f1[i],f2[i],fwhm[i]))
    f.close()
    print("File new_sizes0.dat written to disk.")



