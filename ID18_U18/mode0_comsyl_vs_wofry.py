import numpy
from srxraylib.plot.gol import plot

c_h = numpy.loadtxt("data/comsyl7_h.txt")
c_v = numpy.loadtxt("data/comsyl7_v.txt")
w_h = numpy.loadtxt( "data/wofry7_h.txt")
w_v = numpy.loadtxt( "data/wofry7_v.txt")
g_h = numpy.loadtxt(   "data/gsm7_h.txt")
g_v = numpy.loadtxt(   "data/gsm7_v.txt")

plot(c_h[:,0], c_h[:,1] / c_h[:,1].max(),
     w_h[:,0], w_h[:,1] / w_h[:,1].max(),
     g_h[:,0], g_h[:,1] / g_h[:,1].max(),
     legend=["COMSYL H", "WOFRY H", "GSM H"])

plot(c_v[:,0], c_v[:,1] / c_v[:,1].max(),
     w_v[:,0], w_v[:,1] / w_v[:,1].max(),
     g_v[:,0], g_v[:,1] / g_v[:,1].max(),
     legend=["COMSYL V", "WOFRY V","GSM V"])