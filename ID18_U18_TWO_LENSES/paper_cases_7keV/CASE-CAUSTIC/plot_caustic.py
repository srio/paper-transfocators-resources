import numpy
from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show
import matplotlib.pylab as plt
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import h5py
from oasys.util.oasys_util import get_fwhm
from barc4plots.barc4plots import Image2Plot, ESRF_colors_2D

fontsize = 21
# fontsize_legend = 22
# matplotlib.rc('xtick', labelsize=fontsize)
# matplotlib.rc('ytick', labelsize=fontsize)
# params = {'legend.fontsize':     fontsize,
#           'legend.handlelength': fontsize // 20}
# plt.rcParams.update(params)


plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'figure.autolayout': True})

plt.rc('font', size=fontsize)
# plt.rc('axes', titlesize=fontsize)
plt.rc('axes', labelsize=fontsize * 8 // 10)
plt.rc('xtick', labelsize=fontsize * 8 // 10)
plt.rc('ytick', labelsize=fontsize * 8 // 10)


create_h5 = 0
directories = [
                'case1h', # 'case1v',
                'case2h', # 'case2v',
                'case3h', # 'case3v',
                'case4h', # 'case4v',
                'case1v',
                'case2v',
                'case3v',
                'case4v',
               ]

# directories = ['case2h0', 'case4v0', 'case4vCorrector', 'case2hCorrector']

# directories = ['case2hRefined', 'case4vRefined',]

F1 = [46.1, 25.1, 46.1, 25.1,
      15.0, 42.2, 85.2, 42.2]
F1size = []



for dd,directory in enumerate(directories):
    if create_h5:
        distances = numpy.linspace(10.0, 50.0, 81)
        for i, distance in enumerate(distances):
            a = numpy.loadtxt("%s/%s_wofry_spectral_density_%g.dat" % (directory, directory, distance))
            print(a.shape)
            if i == 0:
                mesh = numpy.zeros(( distances.size, a.shape[0] ))
            mesh[i,:] = a[:,1]


        wr = H5SimpleWriter.initialize_file(filename=directory+".h5", creator="srio", overwrite=1)
        wr.create_entry("caustic", nx_default="Intensity")
        wr.add_image(mesh, distances, a[:,0], entry_name="caustic",image_name="Intensity",
                            title_x="distance [m]",title_y=r'X [$\mu$m]')


    f = h5py.File("%s.h5" % directory,'r')
    mesh = f["caustic/Intensity/image_data"][()].T
    distances = f["caustic/Intensity/axis_x"][()]
    y = f["caustic/Intensity/axis_y"][()]
    f.close()





    mesh_central = mesh[:, mesh.shape[1] // 2]
    ii = numpy.argmax(mesh_central)

    FWHM = numpy.zeros_like(distances)
    for i in range(distances.size):
        fwhm, quote, coordinates = get_fwhm(mesh[i, :], y)
        FWHM[i] = fwhm


    do_plot = 0

    if do_plot:
        # plot_image(mesh, distances - 30.0, y, aspect='auto',
        #            title=directory, xtitle="distance from focus [m]", ytitle="spatial coordinate [um]", show=0)
        #
        #
        # plot(y, mesh[mesh.shape[0]//2, :],
        #      y, mesh[ii, :], legend=['central', 'best focus'], show=0)
        #
        # plot(distances - 30.0, mesh_central, show=0)
        #
        # plot(distances - 30.0, FWHM, xtitle="distance from focus [m]", ytitle="FWHM [um]")

        plot_image_with_histograms(mesh, distances - 30.0, y, aspect='auto', figsize=(10,10),
                title=directory, xtitle="distance from focus [m]", ytitle=r'spatial coordinate [$\mu$m]',
                use_profiles_instead_histograms=1, show=0,
                cmap=ESRF_colors_2D(8),)

        filepng = "%s_caustic.png" % directory
        plt.savefig(filepng)
        print("File written to disk: %s" % filepng)

        plot_show()


    ibest = numpy.argmin(FWHM)
    fwhm_limit = FWHM[ibest] * 1.25
    i_limit = numpy.argwhere(FWHM < fwhm_limit)
    # print(">>>>", i_limit, distances[i_limit[0]], distances[i_limit[-1]])
    dof = distances[i_limit[-1]] - distances[i_limit[0]]
    print("%s best focus: %g um at %g m (depth: %g)" % (directory, FWHM[ibest], distances[ibest] - 30, dof))
    # print(">>>> %g %g" % (F1[dd], FWHM[ibest]) )
    F1size.append((FWHM[ibest]))



print(F1)
print(F1size)