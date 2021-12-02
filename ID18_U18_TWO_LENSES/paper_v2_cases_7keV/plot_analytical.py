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

    M = (q1 / q1) * (q2 / p2)
    return f2, M

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

    fig2, ax2 = plot(
         numpy.array(F1), Msource_at_id * sourcesize_h * 1e6,
         numpy.array(F1), Msource_at_slit * aperture_h * 1e6,
         numpy.array(F1), Msource_at_id * sourcesize_v * 1e6,
         numpy.array(F1), Msource_at_slit * aperture_v * 1e6,
         color=['b','b','r','r'],
         marker=[None,None,None,None],
         linestyle=['--',":",'--',":"],
         yrange=[-5,1e2],
         ylog=0,
         legend=["H (Source at ID)","H (Source at slit)","V (Source at ID)","V (Source at slit)"],
         xtitle="F1 [m]", ytitle="FWHM [$\mu$m]", title="Sizes. Slit=H x V= %g x %g $\mu$m$^2$" %
                                                    (aperture_h*1e6, aperture_v*1e6),
         show=0)

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

    plot_sizes_analytical(aperture_h=40.3e-6, aperture_v=227e-6, filename='_analytical.eps')


