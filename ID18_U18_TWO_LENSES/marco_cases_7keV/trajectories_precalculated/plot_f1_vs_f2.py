import numpy
from srxraylib.plot.gol import plot

a_2h = numpy.loadtxt("f1_vs_f2_case2h.dat")
a_2v = numpy.loadtxt("f1_vs_f2_case2v.dat")

a_8h = numpy.loadtxt("f1_vs_f2_case8h.dat")
a_8v = numpy.loadtxt("f1_vs_f2_case8v.dat")

plot(a_2h[:,0], a_2h[:,1],
     a_2v[:,0], a_2v[:,1],
     a_8h[:,0], a_8h[:,1],
     a_8v[:,0], a_8v[:,1],
     legend=["case2h", "case2v", "case8h", "case8v",],
     title="trajectories")