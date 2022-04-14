import numpy
from srxraylib.plot.gol import plot, plot_image
from silx.io.specfile import SpecFile
from oasys.util.oasys_util import get_fwhm
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import h5py
from plot_analytical import get_f2

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


    flag_create_dat_file = False # use False for creating images, True for creating dat file


    direction = 'H'

    #
    #  theoretical
    #

    if flag_create_dat_file:
        F1 = numpy.linspace(5, 100, 96)
    else:
        F1 = numpy.linspace(5, 100, 500)




    # F2 = numpy.linspace(1, 61, 183)


    F2_Source = numpy.zeros_like(F1)
    F2_Slit = numpy.zeros_like(F1)


    for index, f1 in enumerate(F1):
        ff_id, mm_source_at_id = get_f2(f1=f1,
                                        position_source=0.0,  # source at source
                                        position_lens1=65.0,
                                        position_lens2=170.0,
                                        position_sample=200.0,
                                        verbose=False)
        ff_slit, mm_source_at_slit = get_f2(f1=F1[index],
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

        F2_Source[index] = ff_id
        F2_Slit[index] = ff_slit

    #
    # matplotlib settings
    #

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



    if direction == 'V':

        APERTURE = [25e-6, 227.0e-6, 506.7e-6, 1500e-6]   # vertical

        for i,aperture in enumerate(APERTURE):

            if flag_create_dat_file:
                filedat = "plot_combined_results_%s_%d.dat" % (direction, i)
                ff = open(filedat, 'w')
                ff.write("# aperture           f1       f2_source   f2_slit   f2_numeric   fwhm\n")

            d0 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_v.h5", aperture)
            # d1 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_firstbranch_v.h5", aperture)
            # d2 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_secondbranch_v.h5", aperture)

            fig, ax = plot_image(d0["F1F2_I0"], d0["F1"], d0["F2"],
                       title="$a_v$=%g $\mu$m" % (1e6*aperture),
                       xtitle="$f_1$ [m]", ytitle="$f_2$ [m]",
                       add_colorbar=0, figsize=(8,8),
                       aspect='auto', show=0)

            IBEST = numpy.zeros_like(d0["F1"], dtype=int)
            for j in range(IBEST.size):
                IBEST[j] = numpy.argmax(d0["F1F2_I0"][j,:])
            # print(">>>",d0["ibest"], IBEST,d0["ibest"].shape, IBEST.shape)

            ax.plot(F1, F2_Source, color='white', linestyle=(0, (1, 10)))
            ax.plot(F1, F2_Slit, color='white', linestyle=(0, (5, 10)))
            # ax.plot(d0["F1"][0::4], d0["F2"][d0["ibest"]][0::4], color='red', linestyle="", marker='.')
            ax.plot(d0["F1"][0::4], d0["F2"][IBEST][0::4], color='red', linestyle="", marker='.')


            filename = "V_%d.png" % i
            plt.savefig(filename)
            print("File written to disk: %s" % filename)

            # plot(d0["F1"], d0["F2"][d0["ibest"]],
            #      # d1["F1"], d1["F2"][d1["ibest"]],
            #      # d2["F1"], d2["F2"][d2["ibest"]],
            #      title="trajectories", xtitle="f1", ytitle="f2", show=0)
            #
            # plot(d0["F1"], d0["fwhm_best"],
            #      # d1["F1"], d1["fwhm_best"],
            #      # d2["F1"], d2["fwhm_best"],
            #      title="beam size", xtitle="f1", ytitle="fwhm", ylog=1, show=0) #, yrange=[1e0,1e3])

            plt.show()

            if flag_create_dat_file:
                for iii in range(F1.size):
                    ff.write("%11g  %11g  %11g %11g  %11g  %11g\n" % (1e6*APERTURE[i], F1[iii], F2_Source[iii], F2_Slit[iii],
                                                                d0["F2"][IBEST][iii], d0["F1F2_FWHM"][iii,IBEST[iii]] ))
                ff.close()
                print("File written to disk: ", filedat)


    elif direction == 'H':
        APERTURE = [40.3e-6, 85.1e-6, 145.5e-6, 1000e-6]  # horizontal

        for i,aperture in enumerate(APERTURE):

            if flag_create_dat_file:
                filedat = "plot_combined_results_%s_%d.dat" % (direction, i)
                ff = open(filedat, 'w')
                ff.write("# aperture           f1       f2_source   f2_slit   f2_numeric   fwhm\n")


            d0 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_h.h5", aperture)
            # d0 = get_data_from_h5("/scisoft/users/srio/COMSYL-SLURM/trajectories_precalculated_mapping/f1f2map_h.h5", aperture)
            # d1 = get_data_from_h5("/scisoft/users/srio/data_id18/trajectories_precalculated_mapping/f1f2map_h_old.h5", aperture)

            print(">>>>>>>>>>>>>>>>>>>>>>",numpy.nanmin(d0["F1F2_I0"]))
            d0["F1F2_I0"] = numpy.nan_to_num(d0["F1F2_I0"], nan=numpy.nanmin(d0["F1F2_I0"]))

            for key in d0.keys():
                print(">>>>", key)

            fig , ax = plot_image(d0["F1F2_I0"], d0["F1"], d0["F2"],
                       title="$a_h$=%g $\mu$m" % (1e6*aperture),
                       xtitle="$f_1$ [m]", ytitle="$f_2$ [m]",
                       add_colorbar=0, figsize=(8, 8),
                       aspect='auto', show=0)


            IBEST = numpy.zeros_like(d0["F1"], dtype=int)
            for j in range(IBEST.size):
                IBEST[j] = numpy.argmax(d0["F1F2_I0"][j,:])
            # print(">>>",d0["ibest"], IBEST,d0["ibest"].shape, IBEST.shape)

            ax.plot(F1, F2_Source, color='white', linestyle=(0, (1, 10)))
            ax.plot(F1, F2_Slit, color='white', linestyle=(0, (5, 10)))
            # ax.plot(d0["F1"][0::4], d0["F2"][d0["ibest"]][0::4], color='red', linestyle="", marker='.')
            ax.plot(d0["F1"][0::4], d0["F2"][IBEST][0::4], color='red', linestyle="", marker='.')

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

            if flag_create_dat_file:
                for iii in range(F1.size):
                    ff.write("%11g  %11g  %11g %11g  %11g  %11g\n" % (1e6*APERTURE[i], F1[iii], F2_Source[iii], F2_Slit[iii],
                                                                d0["F2"][IBEST][iii], d0["F1F2_FWHM"][iii,IBEST[iii]] ))
                ff.close()
                print("File written to disk: ", filedat)

