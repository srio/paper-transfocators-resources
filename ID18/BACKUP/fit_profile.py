import numpy
from srxraylib.plot.gol import plot
a = numpy.loadtxt("/users/srio/Oasys/myprofile2.txt")
print(a.shape)

n = a.shape[0]

w =  n // 20


plot(a[:,0], a[:,1],
     a[(n // 2 - w):(n // 2 + w),0], a[(n // 2 - w):(n // 2 + w),1])

x = a[(n // 2 - w):(n // 2 + w),0] * 1e-6
y = a[(n // 2 - w):(n // 2 + w),1] * 1e-6

yder = numpy.gradient(y, x)
coeff = numpy.polyfit(x, yder, 1)
plot(x,yder, x, coeff[0]*x + coeff[1])
print("lens (with two curved sides) of  radius = %g um " % (1e6 * 2/coeff[0]))

