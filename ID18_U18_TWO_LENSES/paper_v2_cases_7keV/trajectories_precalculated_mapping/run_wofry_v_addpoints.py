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

                if numpy.isnan(ff_id) or numpy.isnan(mm_source_at_id) or numpy.isnan(ff_slit) or numpy.isnan(mm_source_at_slit):
                    pass # do not calculate this case
                else:
                    expected_size = sourcesize_v * (1-q2/F2[j]) / (1-p1/F1[i])
                    size1 = sourcesize_v * mm_source_at_id
                    size2 = aperture * mm_source_at_slit
                    size_min = numpy.min([size1, size2])
                    size_max = numpy.max([size1, size2])
                    interval = 0.25 * numpy.abs(size1 - size2)

                    size_mean = 0.5 * (size1 + size2)
                    if expected_size > (size_min - interval) and expected_size < (size_max + interval):
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
