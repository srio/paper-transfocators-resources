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

def plot_sizes_analytical(aperture_h=10e-6, aperture_v=10e-6, filename=None):

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


    print("F:", len(F1))


    F2theory1smooth = numpy.array(F2theory1)
    F2theory2smooth = numpy.array(F2theory2)
    Msource_at_id   = numpy.abs(numpy.array(Msource_at_id))
    Msource_at_slit = numpy.abs(numpy.array(Msource_at_slit))

###############################################
    if True:
        fig, ax = plot(
                       F1, numpy.array(ff_source_at_id),
                       F1, numpy.array(ff_source_at_slit),
                       color=['k','k'],
                       linestyle=['--',':'],
                       marker=[None, None],
                       legend=['Source at ID', 'Source at slit'],
                       yrange=[0, 60],
                       # title="trajectories",
                       xtitle=r'$f_1$ [m]', ytitle=r'$f_2$ [m]', title="",
                       show=0)

        ax.xaxis.grid()
        ax.yaxis.grid()

        if filename is not None:
            plt.savefig('f1f2'+filename)
            print("File written to disk: f1f2%s" % filename)

        plt.show()

################################################
    if True:
        fig2, ax2 = plot(numpy.array(F1), numpy.abs(Msource_at_id),
             numpy.array(F1), numpy.abs(Msource_at_slit),
             color=['k','k'], linestyle=['--',':'],
             legend=['Source at ID', 'Source at slit'],
             xtitle="F1 [m]", ytitle="M", ylog=1, yrange=[1e-3,1e2],
             show=0)


        ax2.xaxis.grid()
        ax2.yaxis.grid()



        if filename is not None:
            plt.savefig('magnification'+filename)
            print("File written to disk: magnification%s" % filename)

        plt.show()


    ############################
    # sizes

    if True:
        Hsrc = Msource_at_id * sourcesize_h * 1e6
        Hslt = Msource_at_slit * aperture_h * 1e6
        Vsrc = Msource_at_id * sourcesize_v * 1e6
        Vslt = Msource_at_slit * aperture_v * 1e6
        fig2, ax2 = plot(
             numpy.array(F1), Hsrc,
             numpy.array(F1), Hslt,
             numpy.array(F1), Vsrc,
             numpy.array(F1), Vslt,
             # numpy.array(F1), Hsrc *          (0.13 / 0.7) + Hslt * (1 - 0.13 / 0.7),
             # numpy.array(F1), Vsrc * numpy.abs(0.58 / 0.7) + Vslt * (1 - 0.58 / 0.7),
             # numpy.array(F1), Hsrc *          (0.13 / 0.9) + Hslt * (1 - 0.13 / 0.9),
             # numpy.array(F1), Vsrc * numpy.abs(0.58 / 0.9) + Vslt * (1 - 0.58 / 0.9),
             color=['b','b','r','r','cyan','pink','cyan','pink'],
             marker=[None,None,None,None,None,None,None,None],
             linestyle=['--',":",'--',":", ':', ':', None, None],
             yrange=[1e0,1e3], #[0,100], #[-5,1e2],
             ylog=1,
             legend=["H (Source at ID)","H (Source at slit)","V (Source at ID)","V (Source at slit)",
                     "H scaled CF=0.7","V scaled CF=0.7", "H scaled CF=0.9","V scaled CF-0.9"],
             xtitle="F1 [m]", ytitle="FWHM [$\mu$m]", title="Sizes. Slit=H x V= %g x %g $\mu$m$^2$" %
                                                        (aperture_h*1e6, aperture_v*1e6),
             show=0)

        ax2.xaxis.grid()
        ax2.yaxis.grid()


        if filename is not None:
            plt.savefig('sizes'+filename)
            print("File written to disk: sizes%s" % filename)

        plt.show()




#
#
#
def plot_sizes_analytical_and_numerical(filename=None):


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



    if False:

        # if direction == 'h':
        #     APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]
        # else:
        #     APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]

        # HH = numpy.loadtxt("plot_combined_results_H_0.dat", skiprows=1)
        # VV = numpy.loadtxt("plot_combined_results_V_0.dat", skiprows=1)
        # CF = [0.9, 0.999]
        # aperture_h = 40.3e-6
        # aperture_v = 25e-6


        # HH = numpy.loadtxt("plot_combined_results_H_0.dat", skiprows=1)
        # VV = numpy.loadtxt("plot_combined_results_V_1.dat", skiprows=1)
        # CF = [0.9, 0.9]
        # aperture_h = 40.3e-6
        # aperture_v = 227.0e-6
        #
        #
        HH = numpy.loadtxt("plot_combined_results_H_1.dat", skiprows=1)
        VV = numpy.loadtxt("plot_combined_results_V_2.dat", skiprows=1)
        CF = [0.7, 0.7]
        aperture_h = 85.1e-6
        aperture_v = 506.7e-6
        # #
        #
        # HH = numpy.loadtxt("plot_combined_results_H_3.dat", skiprows=1)
        # VV = numpy.loadtxt("plot_combined_results_V_3.dat", skiprows=1)
        # CF = [0.13, 0.58]
        # aperture_h = 1000e-6
        # aperture_v = 1500e-6


        Hsrc = Msource_at_id * sourcesize_h * 1e6
        Hslt = Msource_at_slit * aperture_h * 1e6
        Vsrc = Msource_at_id * sourcesize_v * 1e6
        Vslt = Msource_at_slit * aperture_v * 1e6

        # linear interpolation
        m_H = (Hslt - Hsrc) / (1 - 0.13)
        m_V = (Vslt - Vsrc) / (1 - 0.58)
        y_H = m_H * (CF[0] - 0.13) + Hsrc
        y_V = m_V * (CF[1] - 0.58) + Vsrc


        fig2, ax2 = plot(
             numpy.array(F1), Hsrc,
             numpy.array(F1), Hslt,
             numpy.array(F1), Vsrc,
             numpy.array(F1), Vslt,
             numpy.array(F1),  y_H, # Hsrc *          (0.13 / CF[0]) + Hslt * (1 - 0.13 / CF[0]),
             numpy.array(F1),  y_V, # Vsrc * numpy.abs(0.58 / CF[1]) + Vslt * (1 - 0.58 / CF[1]),
             # numpy.array(F1), Hsrc *          (0.13 / 0.9) + Hslt * (1 - 0.13 / 0.9),
             # numpy.array(F1), Vsrc * numpy.abs(0.58 / 0.9) + Vslt * (1 - 0.58 / 0.9),
             HH[:, 1], HH[:, -1],
             VV[:, 1], VV[:, -1],
             color=['b','b','r','r','cyan','pink','cyan','pink'],
             marker=[None,None,None,None,None,None,'.','.'],
             linestyle=['--',":",'--',":", None, None,'',''],
             yrange=[1,1e3], #[0,100], #[-5,1e2],
             ylog=1,
             legend=["H (Source at ID)","H (Source at slit)","V (Source at ID)","V (Source at slit)",
                     "H scaled CF=%g" % CF[0], "V scaled CF=%g" % CF[1],
                     "numeric H", "numeric V"],
             xtitle="F1 [m]", ytitle="FWHM [$\mu$m]", title="Sizes. Slit=H x V= %g x %g $\mu$m$^2$" %
                                                        (1e6*aperture_h, 1e6*aperture_v),
             show=0)

        # ax2.fill_between(numpy.array(F1), Hsrc, Hslt,)
        # ax2.fill_between(numpy.array(F1), Vsrc, Vslt, )
        ax2.xaxis.grid()
        ax2.yaxis.grid()

        if filename is not None:
            plt.savefig('sizes'+filename)
            print("File written to disk: sizes%s" % filename)

        plt.show()

    if True:


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
        fig2, ax2 = plot(
            HH0[:, 1], HH0[:, -1],
            HH1[:, 1], HH1[:, -1],
            HH2[:, 1], HH2[:, -1],
            HH3[:, 1], HH3[:, -1],
            color=['b','r','cyan','pink',],
             marker=['.','.','.','.',],
             linestyle=['','','',''],
             # yrange=[1,1e3], ylog=1,
             yrange=[0,50], ylog=0,
             legend=["a=%g" % (1e6*APERTURE_H[0]),"a=%g" % (1e6*APERTURE_H[1]),"a=%g" % (1e6*APERTURE_H[2]),"a=%g" % (1e6*APERTURE_H[3])],
             xtitle="F1 [m]", ytitle="FWHM [$\mu$m]", title="Sizes H",
             show=0)

        ax2.fill_between(HH0[:, 1], HH0[:, -1], HH3[:, -1])
        ax2.xaxis.grid()
        ax2.yaxis.grid()

        if filename is not None:
            plt.savefig('sizes'+filename)
            print("File written to disk: sizes%s" % filename)

        plt.show()

        # v
        fig2, ax2 = plot(
            VV0[:, 1], VV0[:, -1],
            VV1[:, 1], VV1[:, -1],
            VV2[:, 1], VV2[:, -1],
            VV3[:, 1], VV3[:, -1],
            color=['b','r','cyan','pink'],
             marker=['.','.','.','.'],
             linestyle=['','','','',],
             # yrange=[1,1e3], ylog=1,
             yrange=[0,50], ylog=0,
             legend=["a=%g" % (1e6*APERTURE_V[0]),"a=%g" % (1e6*APERTURE_V[1]),"a=%g" % (1e6*APERTURE_V[2]),"a=%g" % (1e6*APERTURE_V[3])],
             xtitle="F1 [m]", ytitle="FWHM [$\mu$m]", title="Sizes V",
             show=0)

        ax2.fill_between(VV0[:, 1], VV0[:, -1], VV3[:, -1])
        ax2.xaxis.grid()
        ax2.yaxis.grid()

        if filename is not None:
            plt.savefig('sizes'+filename)
            print("File written to disk: sizes%s" % filename)

        plt.show()
if __name__ == "__main__":


    # if direction == 'h':
    #     APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]
    # else:
    #     APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]

    # plot_sizes_analytical(aperture_h=40.3e-6, aperture_v=227e-6, filename='_analytical.eps')
    # plot_sizes_analytical(aperture_h=1000e-6, aperture_v=1500e-6, filename='_analytical.eps')


    plot_sizes_analytical_and_numerical()

