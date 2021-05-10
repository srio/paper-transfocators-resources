from srxraylib.plot.gol import plot, set_qt
import matplotlib.pylab as plt

set_qt()
#
import h5py
import numpy
FACTOR = [0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5] #, 2, 4, 6]
FACTOR = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5] #, 2, 4, 6]


MINFWHM_UNDSOURCE_RECTSLIT_R200UM = []
MINFWHM_GAUSSSOURCE_RECTSLIT_R200UM = []
MINFWHM_GAUSSSOURCE_GAUSSSLIT_R200UM = []
MINFWHM_GAUSSSOURCE_GAUSSSLIT_IDEALLENS = []
MINFWHM_HYBRID = []

for i in range(len(FACTOR)):
        filein = "%s/aperture_factor_%g.dat" % ("./7keV_UndSource_RectSlit_R200um", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM_UNDSOURCE_RECTSLIT_R200UM.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])

        #

        filein = "%s/aperture_factor_%g.dat" % ("./7keV_GaussianSource_RectSlit_R200um", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM_GAUSSSOURCE_RECTSLIT_R200UM.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])

        #

        filein = "%s/aperture_factor_%g.dat" % ("./7keV_GaussianSource_GaussianSlit_R200um", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM_GAUSSSOURCE_GAUSSSLIT_R200UM.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])




        filein = "%s/aperture_factor_%g.dat" % ("./7keV_GaussianSource_GaussianSlit_IdealLens", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM_GAUSSSOURCE_GAUSSSLIT_IDEALLENS.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])


        #
        filein = "%s/rt_aperture_factor_%g.h5" % ("./7keV_ShadowHybrid_R200um", FACTOR[i])
        a1 = h5py.File(filein, 'r')
        # print(a1.shape)
        distance1 = numpy.array(a1['/caustic/fwhm/x'])
        fwhm1 = numpy.array(a1['/caustic/fwhm/y'])
        icenter1 = numpy.array(a1['/caustic/I0/y'])
        MINFWHM_HYBRID.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])



beam_dimension_at_slit_in_um = 565

plot(numpy.array(FACTOR),numpy.array(MINFWHM_UNDSOURCE_RECTSLIT_R200UM),
     numpy.array(FACTOR), numpy.array(MINFWHM_GAUSSSOURCE_RECTSLIT_R200UM),
     numpy.array(FACTOR), numpy.array(MINFWHM_GAUSSSOURCE_GAUSSSLIT_R200UM),
     numpy.array(FACTOR), numpy.array(MINFWHM_GAUSSSOURCE_GAUSSSLIT_IDEALLENS),
     numpy.array(FACTOR),numpy.array(MINFWHM_HYBRID),
     numpy.array(FACTOR)[-2:],numpy.array(MINFWHM_HYBRID)[-2:] * 0 + 9.13 * (18.8 / 65),
     numpy.array(FACTOR)[0:2],( numpy.array(FACTOR)[0:2] * beam_dimension_at_slit_in_um) * (28.8 / (65-35) ),
     xtitle="n", ytitle="WAIST FWHM [um]", ylog=0,
     marker=["o","o","o","o","o", "None","None"],
     legend=["U source + R slit + R=200um",
             "G* source + R slit + R=200um",
             "G* source + G slit + R=200um",
             "G* source + G slit + IdealL",
             "Hybrid",
             "GeomOpt at source",
             "GeomOpt at slit"],
     show=0)

file_png = "waist_size_comparison.eps"
plt.savefig(file_png)
print("File written to disk: %s" % file_png)

plt.show()