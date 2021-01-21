

import numpy
from srxraylib.plot.gol import plot


import matplotlib.pylab as plt
def getLandM(f1 = 28.2, f2 = 39.7, verbose=True):



    p1 = position_lens1 - position_source
    q1 = 1 / (1 / f1 - 1 / p1)
    D = position_lens2 - position_lens1
    p2 = D - q1
    q2 = 1.0 / (1 / f2 - 1 / p2)

    L = p1 + q1 + p2 + q2

    M = (1 - q2 / f2) / (1 - p1 / f1)

    if verbose:
        print("D: %g, q1+p2: %g" % (D, q1+p2))
        print("p1: %g" % p1)
        print("q1: %g" % q1)
        print("p2: %g" % p2)
        print("q2: %g" % q2)
        print("L: %g, Sum: %g" % (L, p1+q1+p2+q2))

    return L, M


def get_f2(f1=28.2, L=200, verbose=True):

    p1 = position_lens1 - position_source
    q1 = 1 / (1 / f1 - 1 / p1)
    D = position_lens2 - position_lens1
    p2 = D - q1

    q2 = position_sample - position_lens2
    f2 = 1.0 / (1 / p2 + 1 / q2)

    if verbose:
        print("D: %g, q1+p2: %g" % (D, q1+p2))
        print("p1: %g" % p1)
        print("q1: %g" % q1)
        print("p2: %g" % p2)
        print("q2: %g" % q2)
        print("L: %g, Sum: %g" % (L, p1+q1+p2+q2))

    return f2


def get_f2(f1=28.2, verbose=False):



    p1 = position_lens1 - position_source
    q1 = 1 / (1 / f1 - 1 / p1)
    D = position_lens2 - position_lens1
    p2 = D - q1

    q2 = position_sample - position_lens2
    f2 = 1.0 / (1 / p2 + 1 / q2)

    if verbose:
        print("\n\nD: %g, q1+p2: %g" % (D, q1+p2))
        print("p1: %g" % p1)
        print("q1: %g" % q1)
        print("p2: %g" % p2)
        print("q2: %g" % q2)
        print("Total length: %g" % (p1+q1+p2+q2))

    return f2


if __name__ == "__main__":
    position_source = 0
    position_slit = 35.0
    position_lens1 = 65.0
    position_lens2 = 164.0
    position_sample = 200.0

    #
    # calculate single point
    #
    L, M = getLandM(f1=28.2, f2=39.7)
    print(">>>>> L, M: ", L, M)


    #
    # make plots for paper
    #


    f1 = numpy.linspace(10, 40, 100)

    #
    #
    position_source = 0
    position_slit = 35.0
    position_lens1 = 65.0
    position_lens2 = 164.0
    position_sample = 200.0

    f2a = get_f2(f1=f1, verbose=False)
    La, Ma = getLandM(f1=f1, f2=f2a, verbose=False)


    #
    #
    position_source = 35
    position_slit = 35.0
    position_lens1 = 65.0
    position_lens2 = 164.0
    position_sample = 200.0

    f2aa = get_f2(f1=f1, verbose=False)
    Laa, Maa = getLandM(f1=f1, f2=f2aa, verbose=False)



    #
    #
    position_source = 0
    position_slit = 35.0
    position_lens1 = 65.0
    position_lens2 = 188.0
    position_sample = 200.0

    f2b = get_f2(f1=f1, verbose=False)
    Lb, Mb = getLandM(f1=f1, f2=f2b, verbose=False)

    #
    #
    position_source = 35
    position_slit = 35.0
    position_lens1 = 65.0
    position_lens2 = 188.0
    position_sample = 200.0

    f2bb = get_f2(f1=f1, verbose=False)
    Lbb, Mbb = getLandM(f1=f1, f2=f2bb, verbose=False)


    plot(f1, f2a,
         f1, f2b,
         f1, f2aa,
         f1, f2bb,
         xtitle="f1 [m]", ytitle="f2 [m]", yrange=[0,f2a.max()],
         legend=["D = 99 m", "D = 123 m", "D = 99 m SOURCE AT SLIT", "D = 123 m SOURCE AT SLIT"],
         color=['red','blue','red','blue'], linestyle=[None,None,'--','--'])

    plot(f1, Ma, f1, Mb,
         f1, Maa, f1, Mbb,
         xtitle="f1 [m]", ytitle="magnification", ylog=1, xrange=[f1.min(), 38],
         legend=["D = 99 m", "D = 123 m", "D = 99 m SOURCE AT SLIT", "D = 123 m SOURCE AT SLIT"],
         color=['red','blue','red','blue'], linestyle=[None,None,'--','--'])



    #
    # plot with numerical results
    #
    a = numpy.loadtxt("tmp_uptomode50.dat", skiprows=3)

    plot(f1, f2a,
         f1, f2aa,
         a[:,0], a[:,-1],
         xtitle="f1 [m]", ytitle="f2 [m]", yrange=[0,3*a[:,-1].max()],
         legend=["D = 99 m", "D = 99 m SOURCE AT SLIT", "D = 99 m NUMERIC"],
         color=['red','red','green'], linestyle=[None,'--','-.'])

    plot(f1, Ma,
         f1, Maa,
         a[:,0], a[:,1] / 15.14,
         xtitle="f1 [m]", ytitle="magnification", ylog=1, xrange=[f1.min(), 38],
         legend=["D = 99 m", "D = 99 m SOURCE AT SLIT", "D 99 m NUMERIC"],
         color=['red','red','green'], linestyle=[None,'--','-.'])


    # #
    # # 2D plot
    # #
    # f1 = numpy.linspace(10, 40, 100)
    # f2 = numpy.linspace(10, 40, 200)
    # F1 = numpy.outer(f1, numpy.ones_like(f2))
    # F2 = numpy.outer(numpy.ones_like(f1), f2)
    # L, M = getLandM(f1=F1, f2=F2, verbose=False)
    #
    #
    # L = numpy.clip( L, 100, 300)
    # print("L: ", L, L.max(), L.min())
    #
    #
    # from srxraylib.plot.gol import plot_surface, plot_image, plot_contour
    #
    # fig = plot_contour(L, f1, f2, title="L", fill=True, contour_levels=100, show=0)
    #
    # plt.contour(f1, f2, L.T,
    #                   levels=[200],
    #                   #colors=['#808080', '#A0A0A0', '#C0C0C0'],
    #                   extend='both')
    #
    # plt.show()




