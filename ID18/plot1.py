import numpy
from srxraylib.plot.gol import plot
a1 = numpy.loadtxt("tmp_h_mode0.dat")
DISTANCE1 = a1[:,0]
FWHM1 = a1[:,1]
I01 = a1[:,2]
ITOTAL1 = a1[:,3]

a2 = numpy.loadtxt("tmp_h_50modes.dat")
DISTANCE2 = a2[:,0]
FWHM2 = a2[:,1]
I02 = a2[:,2]
ITOTAL2 = a2[:,3]

a3 = numpy.loadtxt("tmp_v_mode0.dat")
DISTANCE3 = a3[:,0]
FWHM3 = a3[:,1]
I03 = a3[:,2]
ITOTAL3 = a3[:,3]

a4 = numpy.loadtxt("tmp_v_5modes.dat")
DISTANCE4 = a4[:,0]
FWHM4 = a4[:,1]
I04 = a4[:,2]
ITOTAL4 = a4[:,3]

# plot(DISTANCE, I0, ytitle="Peak intensity [a.u.]", xtitle="Distance from source [m]", figsize=(15, 4), show=0)
plot(
        DISTANCE1, 1e6*FWHM1,
        DISTANCE2, 1e6*FWHM2,
        DISTANCE3, 1e6 * FWHM3,
        DISTANCE4, 1e6 * FWHM4,
        ytitle="FWHM [um]", xtitle="Distance from source [m]", figsize=(15, 4), show=0,
        legend=["H mode 0", "H 50 modes","V mode 0", "V 5 modes"],
        linestyle=['--',None,'--',None],
        color=['red','red','blue','blue'])

plot(
        DISTANCE1, I01,
        DISTANCE2, I02,
        DISTANCE3, I03,
        DISTANCE4, I04,
        ytitle="Intensity at center", xtitle="Distance from source [m]", figsize=(15, 4), show=0,
        legend=["H mode 0", "H 50 modes", "V mode 0", "V 5 modes"],
        linestyle=['--', None, '--', None],
        color=['red','red','blue','blue'])

plot(
        DISTANCE1, ITOTAL1,
        DISTANCE2, ITOTAL2,
        DISTANCE3, ITOTAL3,
        DISTANCE4, ITOTAL4,
        ytitle="Beam Intensity", xtitle="Distance from source [m]", figsize=(15, 4), show=1,
        legend=["H mode 0", "H 50 modes", "V mode 0", "V 5 modes"],
        linestyle=['--', None, '--', None],
        color=['red','red','blue','blue'])


FWHM1[0:(FWHM1.size//2)] = 1
FWHM2[0:(FWHM2.size//2)] = 1
FWHM3[0:(FWHM3.size//2)] = 1
FWHM4[0:(FWHM4.size//2)] = 1

iMin1 = numpy.argmin(FWHM1)
iMin2 = numpy.argmin(FWHM2)
iMin3 = numpy.argmin(FWHM3)
iMin4 = numpy.argmin(FWHM4)

print("Minima found for: H mode 0: %g, H 50 modes: %g, V mode 0: %g, V 5 modes: %g" % (
        DISTANCE1[iMin1],DISTANCE2[iMin2],DISTANCE3[iMin3],DISTANCE4[iMin4] ))
print("Minima values for: H mode 0: %g, H 50 modes: %g, V mode 0: %g, V 5 modes: %g" % (
        1e6*FWHM1[iMin1], 1e6*FWHM2[iMin2], 1e6*FWHM3[iMin3], 1e6*FWHM4[iMin4] ))