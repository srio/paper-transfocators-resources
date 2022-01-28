import numpy
from srxraylib.plot.gol import plot, plot_image
from silx.io.specfile import SpecFile
from oasys.util.oasys_util import get_fwhm
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import h5py

if __name__ == "__main__":

    write_h5 = 0



    APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]  # horizontal
    APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]   # vertical


    # nF1 = 100
    # F1 = numpy.linspace(5, 100, nF1)
    # F2MIN = [15, 15, 10, 5]
    # F2MAX = [35, 35, 50, 55]
    # direction='h'
    # directory='results_h'
    # filename = "f1f2map_%s.h5" % (direction)

    # nF1 = 100
    # F1 = numpy.linspace(5, 100, nF1)
    # F2MIN = [10, 10, 10, 10]
    # F2MAX = [60, 60, 60, 60]
    # direction='v'
    # directory='results'
    # filename = "f1f2map_%s.h5" % (direction)


    # first branch
    nF1 = 41
    F1 = numpy.linspace(5, 45, nF1)
    IN = [10, 15, 10, 10]
    F2MIN = [10, 15, 10, 10]
    F2MAX = [25, 35, 60, 60]
    direction='v'
    directory="results_firstbranch_v"
    filename = "f1f2map_firstbranch%s.h5" % (direction)

    # second branch
    # nF1 = 76
    # F1 = numpy.linspace(25, 100, nF1)
    # F2MIN = [22, 22, 22, 22]
    # F2MAX = [50, 50, 50, 50]
    # direction='v'
    # directory="results_secondbranch_v"
    # filename = "f1f2map_secondbranch_%s.h5" % (direction)


    F2_BEST = numpy.zeros_like(F1)

    if write_h5:
        wr = H5SimpleWriter.initialize_file(filename=filename, creator="srio", overwrite=1)

        for i in range(len(APERTURE)):
            aperture = APERTURE[i]
            F2 = numpy.linspace(F2MIN[i], F2MAX[i], nF1)

            F1F2_FWHM = numpy.zeros((F1.size, F2.size))
            F1F2_I0 = numpy.zeros_like(F1F2_FWHM)
            F1F2_PEAK = numpy.zeros_like(F1F2_FWHM)

            for i in range(F1.size):
                for j,f2 in enumerate(F2):
                    filename="%s/%s_a=%4.2f_f1=%4.2f_f2=%4.2f.dat" % (directory,direction, 1e6*aperture,F1[i],F2[j])

                    try:
                        a = numpy.loadtxt(filename, skiprows=4)

                        fwhm, quote, coordinates = get_fwhm(a[:,1], a[:,0])

                        F1F2_FWHM[i,j] =  float(fwhm)
                        F1F2_I0[i, j] = a[a.shape[0]//2,1]
                        F1F2_PEAK[i, j] = a[:, 1].max()
                    except:
                        pass # print("ERROR " + filename)


            if True:
                plot_image(F1F2_FWHM,F1,F2,title="FWHM aperture=%g" % aperture, aspect='auto')
                plot_image(F1F2_I0,F1,F2,title="I0 aperture=%g" % aperture, aspect='auto')
                plot_image(F1F2_PEAK,F1,F2,title="PEAK aperture=%g" % aperture, aspect='auto')


            entry_name = "f1f2map_a=%4.2f" % (1e6 * aperture)
            wr.create_entry(entry_name, nx_default="I0")

            wr.add_image(F1F2_FWHM,F1,F2, entry_name=entry_name, image_name="FWHM",
                         title_x="f1 [um]", title_y=r'f2 [m]')

            wr.add_image(F1F2_I0,F1,F2, entry_name=entry_name, image_name="I0",
                         title_x="f1 [um]", title_y=r'f2 [m]')

            wr.add_image(F1F2_PEAK,F1,F2, entry_name=entry_name, image_name="PEAK",
                         title_x="f1 [um]", title_y=r'f2 [m]')


            wr.create_entry("%s/scans" % entry_name) #, nx_default="aaa")

            for i in range(F1.size):
                for j,f2 in enumerate(F2):
                    # sub_entry = "%s/scans/f1=%4.2f_f2=%4.2f" % (entry_name,F1[i],F2[j] )
                    filename = "results/%s_a=%4.2f_f1=%4.2f_f2=%4.2f.dat" % (direction, 1e6 * aperture, F1[i], F2[j])
                    try:
                        a = numpy.loadtxt(filename, skiprows=4)
                    except:
                        pass # print("ERROR " + filename)

                    wr.add_dataset(a[:,0], a[:,1], "f1=%4.2f_f2=%4.2f" % (F1[i],F2[j] ),
                                   entry_name="%s/scans" % entry_name)

    else:
        f = h5py.File(filename, 'r')
        for aperture in APERTURE:
            entry = "f1f2map_a=%4.2f" % (1e6 * aperture)

            F1F2_FWHM = f["%s/FWHM/image_data" % entry][()].T
            F1F2_I0 = f["%s/I0/image_data" % entry][()].T
            F1F2_PEAK = f["%s/PEAK/image_data" % entry][()].T
            F1 = f["%s/FWHM/axis_x" % entry][()]
            F2 = f["%s/FWHM/axis_y" % entry][()]

            # plot_image(F1F2_FWHM, F1, F2, title="FWHM aperture=%g" % aperture, aspect='auto')
            plot_image(F1F2_I0, F1, F2, title="I0 aperture=%g" % aperture, aspect='auto', show=0)
            # plot_image(F1F2_PEAK, F1, F2, title="PEAK aperture=%g" % aperture, aspect='auto')

            ibest = numpy.argmax(F1F2_I0, axis=1)
            fwhm_best = numpy.zeros_like(F1)
            for i in range(F1.size):
                fwhm_best[i] = F1F2_FWHM[i,ibest[i]]

            plot(F1, F2[ibest], xtitle="f1", ytitle="f2", show=0)
            plot(F1, fwhm_best, xtitle="f1", ytitle="fwhm", ylog=1, yrange=[1e0,1e3])
        f.close()
