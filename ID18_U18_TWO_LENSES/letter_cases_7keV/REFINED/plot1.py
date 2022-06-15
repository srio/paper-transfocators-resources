import numpy
from srxraylib.plot.gol import plot

a = numpy.loadtxt("new_sizes1.dat")
plot(a[:,0], a[:,2],
     a[:,0], a[:,3],
     legend=['old','new'], yrange=[0,50])