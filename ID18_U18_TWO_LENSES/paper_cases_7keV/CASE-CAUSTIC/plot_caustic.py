import numpy
from srxraylib.plot.gol import plot, plot_image
from srxraylib.util.h5_simple_writer import H5SimpleWriter
import h5py

create_h5 = 1
directory = "case4h"


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
           title=directory, xtitle="distance from focus [m]", ytitle="spatial coordinate [um]" )