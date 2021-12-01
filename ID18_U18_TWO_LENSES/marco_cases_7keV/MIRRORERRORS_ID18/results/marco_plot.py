import numpy
import matplotlib.pylab as plt

def go(distance=36, dir="./ID18slpx1/"):
     file_h = '%s/tmp%d_h_spectral_density.dat' % (dir, distance)
     file_v = '%s/tmp%d_v_spectral_density.dat' % (dir, distance)

     h = numpy.loadtxt(file_h)
     v = numpy.loadtxt(file_v)

     h[:,1] /= h[:,1].max()
     v[:,1] /= v[:,1].max()

     hv = numpy.outer(h[:,1], v[:,1])


     z = hv
     x = h[:,0]
     y = v[:,0]

     plt.imshow(z.T, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], cmap=None, aspect='auto')
     plt.title("%s distance=%d"%(dir,distance))
     plt.show()


if __name__ == "__main__":

     try:
          from srxraylib.plot.gol import set_qt
          set_qt()
     except:
          pass

     dir0 = "https://raw.githubusercontent.com/srio/paper-transfocators-resources/main/ID18_U18_TWO_LENSES/marco_cases_7keV/MIRRORERRORS_ID18/results/"

     go(distance=36, dir=dir0 + "/ID18slpx0/")
     go(distance=36, dir=dir0 + "/ID18slpx0p5/")
     go(distance=36, dir=dir0 + "/ID18slpx1/")
     go(distance=36, dir=dir0 + "/ID18slpx2/")