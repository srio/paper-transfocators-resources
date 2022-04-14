import numpy
from datetime import datetime

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

    tally = TallyCoherentModes()
    for my_mode_index in range(10):
        output_wavefront = run_source(my_mode_index=my_mode_index)
        output_wavefront = run_beamline(output_wavefront, aperture, r1, r2)
        tally.append(output_wavefront)

    if do_plot:
        # tally.plot_cross_spectral_density(show=1, filename="")
        tally.plot_spectral_density(show=1, filename="", title="a=%g, r1=%g" % (aperture, r2))
        # tally.plot_occupation(show=1, filename="")

    if filename != "":
        # tally.save_scan(filename)
        # print("Saving file...")
        tally.save_spectral_density(filename)
        # print("File written to disk: %s" % filename)


if __name__ == "__main__":

    import xraylib

    from wofry_beamline_v import run_source, run_beamline

    APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]   # vertical
    sourcesize_v = 15.02e-6
    direction='v'
    directory="results_v"
    p1 = 65.0
    q2 = 30.0

    # write source (the same for all cases)
    for my_mode_index in range(10):
        run_source(my_mode_index=my_mode_index, load_from_file=False, write_to_file=True)

    F1 = numpy.linspace(5, 100, 96)
    F2 = numpy.linspace(1, 61, 183)
    xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
    R1 = F1 * (2 * xrl_delta)
    R2 = F2 * (2 * xrl_delta)

    for ii in range(len(APERTURE)):

        aperture = APERTURE[ii]

        for i in range(F1.size):
            for j in range(F2.size):

                if check_added(F1[i], F2[j], direction, aperture):
                    filename="%s/%s_a=%4.2f_f1=%4.2f_f2=%4.2f.dat" % (directory, direction, 1e6*aperture, F1[i], F2[j])
                    try:
                        f = open(filename, 'r')
                        f.close()
                    except:
                        print("%s i=%d,j=%d,current=%d,total=%d, remains: %g %%" % (
                                                            datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                                                            i,j, i*R1.size+j,
                                                            R1.size*R2.size,100*(1-(i*R1.size+j)/(R1.size*R2.size))) )
                        main(aperture, R1[i], R2[j],
                             filename="%s/%s_a=%4.2f_f1=%4.2f_f2=%4.2f.dat" % (directory, direction, 1e6*aperture, F1[i], F2[j]))
