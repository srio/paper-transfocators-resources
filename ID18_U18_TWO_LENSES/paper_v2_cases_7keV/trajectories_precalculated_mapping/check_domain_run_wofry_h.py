import numpy
from datetime import datetime

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

def main(aperture, r1, r2, filename="", do_plot=0):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes



def check_added(f1, f2, direction, aperture):

    if direction == 'h':
        if aperture == 40.3e-6: # [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]:
            if f1 > 28 and f1 < 32:
                if f2 >= 18.0: return True
            if f1 > 40:
                if f2 >= 21 and f2 <=32: return True
        elif aperture == 85.1e-6: # [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]:
            if f1 == 29.0 and f2 >= 15.0: return True
            if f1 > 40:
                if f2 >= 21 and f2 <=30: return True
        elif aperture == 145.5e-6: # [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]:
            if f1 == 29.0 and f2 >= 15.0: return True
            if f1 > 40:
                if f2 >= 21 and f2 <=30: return True
        elif aperture == 1000e-6: # [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]:
            if f1 == 29.0 and f2 >= 15.0: return True
            if f1 > 57:
                if f2 >= 27 and f2 <=33: return True
    elif direction == 'v':
        if aperture == 25e-6: # [25e-6, 227.0e-6, 506.7e-6, 1500e-6]  # vertical
            if f1 > 28 and f1 < 32:
                if f2 >= 15.0: return True
            if f1 > 40:
                if f2 >= 21 and f2 <=38: return True
        elif aperture == 227.0e-6: # [25e-6, 227.0e-6, 506.7e-6, 1500e-6]  # vertical
            if f1 == 29.0 and f2 >= 15.0: return True
            if f1 > 55:
                if f2 >= 27 and f2 <=33: return True
        elif aperture == 506.7e-6: # [25e-6, 227.0e-6, 506.7e-6, 1500e-6]  # vertical
            if f1 == 29.0 and f2 >= 15.0: return True
            if f1 > 55:
                if f2 >= 27 and f2 <=33: return True
        elif aperture == 1500e-6: # [25e-6, 227.0e-6, 506.7e-6, 1500e-6]  # vertical
            if f1 == 29.0 and f2 >= 15.0: return True
            if f1 > 57:
                if f2 >= 27 and f2 <=33: return True

    return False

if __name__ == "__main__":

    import xraylib
    from srxraylib.plot.gol import plot_image, plot


    direction='v'

    if direction == 'h':
        APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]  # horizontal
        sourcesize = 70.57e-6
    elif direction == 'v':
        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]  # vertical
        sourcesize = 15.02e-6

    p1 = 65.0
    q2 = 30.0


    F1 = numpy.linspace(5, 100, 96)
    F2 = numpy.linspace(1, 61, 183)
    xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
    R1 = F1 * (2 * xrl_delta)
    R2 = F2 * (2 * xrl_delta)

    FF_ID = numpy.zeros_like(F1)
    FF_SLIT = numpy.zeros_like(F1)

    domain = numpy.zeros((len(F1),len(F2)))
    EXPECTED_SIZE = numpy.zeros((len(F1),len(F2)))
    for ii in range(len(APERTURE)):

        aperture = APERTURE[ii]

        for i in range(F1.size):
            for j in range(F2.size):

                ff_id, mm_source_at_id = get_f2(f1=F1[i],
                                                position_source=0.0,  # source at source
                                                position_lens1=65.0,
                                                position_lens2=170.0,
                                                position_sample=200.0,
                                                verbose=False)

                ff_slit, mm_source_at_slit = get_f2(f1=F1[i],
                                                    position_source=36.0,  # 0.0,  # source at slit
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

                FF_ID  [i] = ff_id   # mm_source_at_id   #
                FF_SLIT[i] = ff_slit # mm_source_at_slit #
                if check_added(F1[i], F2[j], direction, aperture):
                    domain[i, j] = 1

                if numpy.isnan(ff_id) or numpy.isnan(mm_source_at_id) or numpy.isnan(ff_slit) or numpy.isnan(mm_source_at_slit):
                # if numpy.isnan(ff_id) or numpy.isnan(mm_source_at_id):
                #     print(">>>>f1: ", F1[i], ff_id, mm_source_at_id)
                    pass # do not calculate this case
                else:
                    expected_size = sourcesize * (1-q2/F2[j]) / (1-p1/F1[i])
                    size1 = sourcesize * mm_source_at_id
                    size2 = aperture * mm_source_at_slit
                    size_min = numpy.min([size1, size2])
                    size_max = numpy.max([size1, size2])
                    interval = 0.25 * numpy.abs(size1 - size2)
                    size_mean = 0.5 * (size1 + size2)
                    if expected_size > (size_min - interval) and expected_size < (size_max + interval):
                        domain[i,j] = 1


        # plot(F1, FF_ID,
        #      F1, FF_SLIT,
        #      yrange=[0, 60], legend=['source at id', 'source at slit'], show=0)
        # plot_image(1e6*EXPECTED_SIZE,F1,F2,title="Aperture %g um" % (1e6 * aperture), show=0)
        plot_image(domain,F1,F2,title="Aperture %g um" % (1e6 * aperture))