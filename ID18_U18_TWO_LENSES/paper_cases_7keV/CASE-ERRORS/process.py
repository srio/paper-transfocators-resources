#
# read files with sizes
#
import numpy
from silx.io.specfile import SpecFile

FWHM = []
X = []
Y = []
# YTOT = []
Index = numpy.arange(200)
case = "1v"
for index in Index:

    filename = "results/case%s_wofry_spectral_density_%03d.dat" % (case, index)
    # filename = "/users/srio/OASYS1.2/paper-transfocators-resources/ID18_U18_TWO_LENSES/paper_cases_7keV/sizes_slit40.3_h/case7keV_h_spectral_density_000.dat"
    print(">>>>>>>>> opening file: ", filename)
    sf = SpecFile(filename)
    s1 = sf[0]
    X.append(s1.data[0, :])
    Y.append(s1.data[1, :])
    # if index == 0:
    #     YTOT = s1.data[1, :]
    # else:
    #     YTOT += s1.data[1, :]
    try:
        fwhm = s1.scan_header_dict["UFWHM"]
        FWHM.append(float(fwhm))
    except:
        FWHM.append(0)

FWHM = numpy.array(FWHM)
print(FWHM.mean(), FWHM.std())

from srxraylib.plot.gol import plot

# plot(
#     X[0], Y[0],
#     X[1], Y[1],
#     X[2], Y[2],
#     X[3], Y[3],
#     X[4], Y[4],
#     X[5], Y[5],
#     X[6], Y[6],
#     X[7], Y[7],
#     X[8], Y[8],
#     X[9], Y[9],
#     X[10], Y[10],
#     X[11], Y[11],
#     X[12], Y[12],
#     title=case
# )

# plot(
#     X[0], YTOT)
