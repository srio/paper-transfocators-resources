import numpy
from srxraylib.plot.gol import plot, plot_image
from silx.io.specfile import SpecFile
from oasys.util.oasys_util import get_fwhm
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import h5py

def get_data_from_h5(filename, aperture):
    f = h5py.File(filename, 'r')
    entry = "f1f2map_a=%4.2f" % (1e6 * aperture)
    F1F2_FWHM = f["%s/FWHM/image_data" % entry][()].T
    F1F2_I0 = f["%s/I0/image_data" % entry][()].T
    F1F2_PEAK = f["%s/PEAK/image_data" % entry][()].T
    F1 = f["%s/FWHM/axis_x" % entry][()]
    F2 = f["%s/FWHM/axis_y" % entry][()]
    f.close()

    ibest = numpy.argmax(F1F2_I0, axis=1)
    fwhm_best = numpy.zeros_like(F1)
    for i in range(F1.size):
        fwhm_best[i] = F1F2_FWHM[i, ibest[i]]


    return {"F1":F1,
            "F2":F2,
            "F1F2_FWHM":F1F2_FWHM,
            "F1F2_I0":F1F2_I0,
            "F1F2_PEAK":F1F2_PEAK,
            "ibest":ibest,
            "fwhm_best":fwhm_best}

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'figure.autolayout': True})

    params = {'legend.fontsize': 15,
              'legend.frameon': False,
              'legend.handlelength': 2,
              'axes.titlesize' : 24,
              'axes.labelsize': 24,
              'lines.linewidth': 3,
              'lines.markersize': 10,
              'xtick.labelsize': 25,
              'ytick.labelsize': 25,
              # 'grid.color':       'r',
              # 'grid.linestyle':   '-',
              # 'grid.linewidth':     2,
              }
    plt.rcParams.update(params)


    direction = 'h'


    if direction == 'v':

        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]   # vertical

        for i,aperture in enumerate(APERTURE):

            d0 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_v.h5", aperture)
            d1 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_firstbranch_v.h5", aperture)
            d2 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_secondbranch_v.h5", aperture)


            plot_image(d0["F1F2_I0"], d0["F1"], d0["F2"],
                       title="$a_v$=%g $\mu$m" % (1e6*aperture),
                       xtitle="$f_1$ [m]", ytitle="$f_2$ [m]",
                       add_colorbar=0, figsize=(8,8),
                       aspect='auto', show=0)


            filename = "V_%d.png" % i
            plt.savefig(filename)
            print("File written to disk: %s" % filename)

            # plot(d0["F1"], d0["F2"][d0["ibest"]],
            #      d1["F1"], d1["F2"][d1["ibest"]],
            #      d2["F1"], d2["F2"][d2["ibest"]],
            #      xtitle="f1", ytitle="f2", show=0)
            #
            # plot(d0["F1"], d0["fwhm_best"],
            #      d1["F1"], d1["fwhm_best"],
            #      d2["F1"], d2["fwhm_best"],
            #      xtitle="f1", ytitle="fwhm", ylog=0, show=0) #, yrange=[1e0,1e3])

            plt.show()

    else:
        APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]  # horizontal

        for i,aperture in enumerate(APERTURE):
            d0 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_h.h5",
                                  aperture)
            d1 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_h_old.h5",
                                  aperture)


            fig , ax = plot_image(d0["F1F2_I0"], d0["F1"], d0["F2"],
                       title="$a_h$=%g $\mu$m" % (1e6*aperture),
                       xtitle="$f_1$ [m]", ytitle="$f_2$ [m]",
                       add_colorbar=0, figsize=(8, 8),
                       aspect='auto', show=0)

            x = numpy.linspace(0,100,100)
            ax.plot(x,x, color='red', linestyle=":")

            filename = "H_%d.png" % i
            plt.savefig(filename)
            print("File written to disk: %s" % filename)

            # plot_image(d1["F1F2_I0"], d1["F1"], d1["F2"],
            #            title="OLD!!!!!!!! %s  I0 aperture=%g" % (direction, aperture),
            #            aspect='auto', show=0)
            #
            # plot(d0["F1"], d0["F2"][d0["ibest"]],
            #      d1["F1"], d1["F2"][d1["ibest"]],
            #      # d2["F1"], d2["F2"][d2["ibest"]],
            #      xtitle="f1", ytitle="f2", show=0)
            #
            # plot(d0["F1"], d0["fwhm_best"],
            #      d1["F1"], d1["fwhm_best"],
            #      # d2["F1"], d2["fwhm_best"],
            #      xtitle="f1", ytitle="fwhm", ylog=0)  # , yrange=[1e0,1e3])

            plt.show()