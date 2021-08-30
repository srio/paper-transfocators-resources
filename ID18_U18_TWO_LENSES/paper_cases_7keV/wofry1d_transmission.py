import numpy
import h5py
from srxraylib.plot.gol import plot_image, set_qt
set_qt()


def get_integral(file, do_plot=1):

    f = h5py.File(file, 'r')

    image = f["/entry/data0/image"][()]
    x = f["/entry/data0/x"][()]
    y = f["/entry/data0/y"][()]

    f.close()


    print(image.shape, x.shape, y.shape)


    if do_plot:
        plot_image(image.T, x, y)


    integral = image.sum() * (x[1] - x[0]) * (y[1] - y[0])



    return integral


if __name__ == "__main__":
    dir = "/Users/srio/Oasys/"
    integral0 = get_integral(dir+"case3_source.h5")
    print(integral0)
    integral1 = get_integral(dir+"case3_slit.h5")
    print(integral1)
    integral2 = get_integral(dir+"case4_lens1.h5")
    print(integral2)
    integral3 = get_integral(dir+"case4_lens2.h5")
    print(integral3)

    print("Absorption slit, lens1, lens2 : %3.1f, %3.1f, %3.1f " % (100 * (integral0 - integral1) / integral0,
                                                      100 * (integral1 - integral2) / integral1,
                                                      100 * (integral2 - integral3) / integral2,) )