import numpy
from srxraylib.plot.gol import plot
import numpy
import h5py
from srxraylib.plot.gol import plot_image, plot, set_qt


def get_wofry_1d_intensities(fileh, filev):
    ah = numpy.loadtxt(fileh)
    av = numpy.loadtxt(filev)

    x = ah[:, 0].copy()
    y = av[:, 0].copy()
    int_x = ah[:, 1].copy()
    int_y = av[:, 1].copy()

    image = numpy.outer(int_x, int_y)

    return image, x, int_x, y, int_y


where = "at_36m"
where = "at_sample"
# dir = "./two_reflections/"
dir = "./one_reflection/"

if True:


    image, x, int_x, y, int_y = get_wofry_1d_intensities(
        dir + "surface_error_%s_h_spectral_density.dat" % where,
        dir + "surface_error_%s_v_spectral_density.dat" % where)


    if False:
        plot(x, int_x,
             y, int_y,
             legend=['h', 'v'])


        plot_image(image, x, y)

    integral = image.sum() * (x[1] - x[0]) * (y[1] - y[0])
    print(integral)

    #load shadow results

    # data saved using plotxy "Save current plot"
    import h5py
    f = h5py.File(dir + "plotxy_%s.hdf5" % where)
    x2 = f["/coordinates/X"][()]
    y2 = f["/coordinates/Y"][()]
    int_x2 = f["/xy_plots/last_plot/histogram_h"][()]
    int_y2 = f["/xy_plots/last_plot/histogram_v"][()]
    f.close()
    print(x2.shape)

    if True:
        plot(x, int_x / int_x.max(),
             1e6 * x2, int_x2 / int_x2.max(),
             legend=['wofry h',
                     'shadow h',
                     ], title=where, xtitle="X [um]", show=0)

        plot(y, int_y / int_y.max(),
             1e6 * y2, int_y2 / int_y2.max(),
             legend=['wofry v',
                     'shadow v',
                     ], title=where, xtitle="Y [um]",)



        # plot_image(image, x, y)




