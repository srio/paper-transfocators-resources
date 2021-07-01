import numpy
from srxraylib.plot.gol import plot
import matplotlib.pylab as plt

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

def plot_f1f2(filename=None,direction='v'):

    if direction == 'h':
        APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]
        LEGEND = ["a=40.3 um", "a=85.1 um", "a=145.5 um", "a=1000 um", "analytical (source at ID)",
                  "analytical (source at slit)"]
    else:
        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]
        LEGEND = ["a=25 um", "a=227 um", "a=506.7 um", "a=1500 um", "analytical (source at ID)",
                  "analytical (source at slit)"]

    a0 = numpy.loadtxt("f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[0], direction), skiprows=3)
    a1 = numpy.loadtxt("f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[1], direction), skiprows=3)
    a2 = numpy.loadtxt("f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[2], direction), skiprows=3)
    a3 = numpy.loadtxt("f1_vs_f2_slit%g_%s.dat" % (1e6 * APERTURE[3], direction), skiprows=3)

    ff_source_at_id, mm_source_at_id = get_f2(f1=a0[:, 0],
                                              position_source=0.0,  # source at source
                                              position_lens1=65.0,
                                              position_lens2=170.0,
                                              position_sample=200.0,
                                              verbose=False)

    ff_source_at_slit, mm_source_at_slit = get_f2(f1=a0[:, 0],
                                                  position_source=36.0,  # 0.0,  # source at slit
                                                  position_lens1=65.0,
                                                  position_lens2=170.0,
                                                  position_sample=200.0,
                                                  verbose=False)

    fig, ax = plot(a0[:, 0], a0[:, 1],
                   a1[:, 0], a1[:, 1],
                   a2[:, 0], a2[:, 1],
                   a3[:, 0], a3[:, 1],
                   a0[:, 0], ff_source_at_id,
                   a0[:, 0], ff_source_at_slit,
                   color=['r', 'b', 'g', 'c', 'k', 'k'],
                   linestyle=[None, None, None, None, '--', ':'],
                   legend=LEGEND,
                   yrange=[0, 60],
                   # title="trajectories",
                   xtitle="f1 [m]", ytitle="f2 [m]",
                   show=0)

    ax.xaxis.grid()
    ax.yaxis.grid()

    if filename is not None:
        plt.savefig(filename)
        print("File written to disk: %s" % filename)

    plt.show()


if __name__ == "__main__":

     plot_f1f2(direction='h', filename="f1f2_h.eps")
     plot_f1f2(direction='v', filename="f1f2_v.eps")

