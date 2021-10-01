import numpy
from srxraylib.plot.gol import plot, plot_image
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import h5py

create_h5 = 0
directories = [
                'case1h', 'case1v',
                'case2h', 'case2v',
                'case3h', 'case3v',
                'case4h', 'case4v',
               ]

create_h5 = 1
directories = ['case2h0', 'case4v0',
               ]

# directory = "case4h"

for directory in directories:
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
                            title_x="distance [m]",title_y="X [um]")


    f = h5py.File("%s.h5" % directory,'r')
    mesh = f["caustic/Intensity/image_data"][()].T
    distances = f["caustic/Intensity/axis_x"][()]
    y = f["caustic/Intensity/axis_y"][()]
    f.close()

    plot_image(mesh, distances - 30.0, y, aspect='auto',
               title=directory, xtitle="distance from focus [m]", ytitle="spatial coordinate [um]", show=0)


    mesh_central = mesh[:, mesh.shape[1]//2]
    ii = numpy.argmax(mesh_central)

    plot(y, mesh[mesh.shape[0]//2, :],
         y, mesh[ii, :], legend=['central', 'best focus'], show=0)

    plot(distances - 30.0, mesh_central)