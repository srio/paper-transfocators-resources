import numpy
from srxraylib.plot.gol import plot, plot_image
from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D

def go(distance=36, do_plot=True,dir="/users/srio/Oasys/ID18slpx1/"):
     file_h = '%s/tmp%d_h_spectral_density.dat' % (dir, distance)
     file_v = '%s/tmp%d_v_spectral_density.dat' % (dir, distance)

     h = numpy.loadtxt(file_h)
     v = numpy.loadtxt(file_v)

     h[:,1] /= h[:,1].max()
     v[:,1] /= v[:,1].max()

     hv = numpy.outer(h[:,1], v[:,1])
     print(hv.shape)

     plot(h[:,0], h[:,1],
          v[:,0], v[:,1],
          title="distance=%d"%distance)

     plot_image(hv,h[:,0],v[:,0],title="distance=%d"%distance)

     wf = GenericWavefront2D.initialize_wavefront_from_arrays(h[:,0]*1e-6, v[:,0]*1e-6,  numpy.sqrt(hv)+0j, z_array_pi=None, wavelength=1e-10)
     wf.save_h5_file("%s/tmp%d.h5" % (dir, distance),subgroupname="wfr",intensity=True,phase=False,overwrite=True,verbose=False)
     print("File witten: ", "%s/tmp%d.h5" % (dir, distance))

go(36)
go(65)
go(200)