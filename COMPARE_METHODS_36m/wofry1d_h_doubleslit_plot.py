from srxraylib.plot.gol import plot, plot_image
import numpy
import h5py




filename = "wofry1d_h_doubleslit.dat"
cf = numpy.loadtxt(filename)
APERTURES = cf[:,0]
CF = cf[:,1]
DoC = cf[:,2]
plot(APERTURES*1e6,CF, xtitle="aperture(outer) [um]", ytitle="CF")
plot(APERTURES*1e6,DoC, xtitle="aperture(outer) [um]", ytitle="DoC fitted")



f = h5py.File("wofry1d_h_doubleslit.h5", 'r')
SD = f["images/Intensity/image_data"][()].T
APERTURES = f["images/Intensity/axis_x"][()]
x = f["images/Intensity/axis_y"][()]
f.close()



plot_image(SD, APERTURES*1e6, x*1e6, xtitle="aperture(outer) [um]",ytitle="abscissas [um]", aspect='auto')

