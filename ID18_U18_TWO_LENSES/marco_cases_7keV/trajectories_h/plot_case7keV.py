import numpy
from srxraylib.plot.gol import plot
from silx.io.specfile import SpecFile

fileroot = "case7keV_2h"
nfiles = 200

#
#
#
Index = numpy.arange(nfiles)
F1 = []
F2 = []
R1 = []
R2 = []
FWHM = []

for index in Index:
    a = numpy.loadtxt("%s_%03d.txt" % (fileroot, index), skiprows=3)
    print(">>>>>> ", a[0], a[1], a[3], a[4])
    R2.append(a[0])
    F2.append(a[1])
    R1.append(a[3])
    F1.append(a[4])

    sf = SpecFile("%s_spectral_density_%03d.dat" % (fileroot, index))
    s1 = sf[0]
    fwhm = s1.scan_header_dict["UFWHM"]
    FWHM.append(float(fwhm))

print("F:", len(F1), len(F2))
print("R:", len(R1), len(R2))

plot(numpy.array(F1), numpy.array(F2), xtitle="F1 [m]", ytitle="F2 [m]", title="trajectories F %s" %fileroot)
plot(1e6 * numpy.array(R1), 1e6 * numpy.array(R2), xtitle="R1 [um]", ytitle="R2 [um]", title="trajectories R %s" %fileroot)
plot(numpy.array(F1), numpy.array(FWHM), xtitle="F1 [m]", ytitle="FWHM [um]", title="Sizes %s" %fileroot)



# filename = "case7keV_2h_spectral_density_000.dat"
#
#
#
#
# # for index in range(len(sf)):
# #     s1 = sf[0]
# name = s1.scan_header_dict["S"]
#
#     # for i,key in enumerate(s1.scan_header_dict.keys()):
#     #     print(i,key,s1.scan_header_dict[key])