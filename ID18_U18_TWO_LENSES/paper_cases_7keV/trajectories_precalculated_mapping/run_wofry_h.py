import numpy
from datetime import datetime

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



    from case2h_wofry import run_source, run_beamline
    :quit()l
    F2MIN = [15, 15, 10, 5]
    F2MAX = [35, 35, 50, 55]
    direction='h'
    directory="results_h"


    # from case2v_wofry import run_source, run_beamline
    # APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]   # vertical
    # F2MIN = [10, 10, 10, 10]
    # F2MAX = [60, 60, 60, 60]
    # direction='v'
    # directory="results_v"


    # for my_mode_index in range(10):
    #     run_source(my_mode_index=my_mode_index, load_from_file=False, write_to_file=True)
    #


    # if True:
    #     ii = 0

    for ii in range(len(APERTURE)):

        aperture = APERTURE[ii]
        f2min = F2MIN[ii]
        f2max = F2MAX[ii]

        F1 = numpy.linspace(5,100,100)
        F2 = numpy.linspace(f2min, f2max, 100)

        xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7, 1.85)).real
        R1 = F1 * (2 * xrl_delta)
        R2 = F2 * (2 * xrl_delta)

        for i in range(R1.size):
            for j in range(R2.size):
                filename="%s/%s_a=%4.2f_f1=%4.2f_f2=%4.2f.dat" % (directory, direction, 1e6*aperture, F1[i], F2[j])
                try:
                    f = open(filename, 'r')
                    f.close()
                except:
                # if True: #F1[i] > 53:
                    print("%s i=%d,j=%d,current=%d,total=%d, remains: %g %%" % (
                                                        datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                                                        i,j, i*R1.size+j,
                                                        R1.size*R2.size,100*(1-(i*R1.size+j)/(R1.size*R2.size))) )
                    main(aperture, R1[i], R2[j],
                         filename="%s/%s_a=%4.2f_f1=%4.2f_f2=%4.2f.dat" % (directory, direction, 1e6*aperture, F1[i], F2[j]))
