from srxraylib.plot.gol import plot, set_qt
set_qt()
#
import h5py
import numpy
FACTOR = [0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5] #, 2, 4, 6]
FACTOR = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5] #, 2, 4, 6]


MINFWHM = []
MINFWHM_HYBRID = []
MINFWHM_GSM = []
MINFWHM_FULLGAUSS = []
MINFWHM_GAUSSSLIT = []
for i in range(len(FACTOR)):
        filein = "%s/aperture_factor_%g.dat" % ("./7keV_R200um_WofryUnd", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])

        #

        filein = "%s/aperture_factor_%g.dat" % ("./7keV_R200um_WofryGSM", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM_GSM.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])

        #

        filein = "%s/aperture_factor_%g.dat" % ("./7keV_R200um_GaussianSlit_WofryGSM", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM_GAUSSSLIT.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])

        #
        filein = "%s/rt_aperture_factor_%g.h5" % ("./7keV_R200um_ShadowHybrid", FACTOR[i])
        a1 = h5py.File(filein, 'r')
        # print(a1.shape)
        distance1 = numpy.array(a1['/caustic/fwhm/x'])
        fwhm1 = numpy.array(a1['/caustic/fwhm/y'])
        icenter1 = numpy.array(a1['/caustic/I0/y'])
        MINFWHM_HYBRID.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])

        filein = "%s/aperture_factor_%g.dat" % ("./7keV_IdealLensGaussianSlit_WofryGSM", FACTOR[i])
        a1 = numpy.loadtxt(filein)
        print(a1.shape)
        distance1 = a1[:,0]  ;  fwhm1 = a1[:,1] ; icenter1 = a1[:,3]
        MINFWHM_FULLGAUSS.append(fwhm1[numpy.argmin(fwhm1)])
        print(">>>>> ", i, filein, fwhm1[numpy.argmin(fwhm1)])

beam_dimension_at_slit_in_um = 565

plot(numpy.array(FACTOR),numpy.array(MINFWHM),
     numpy.array(FACTOR), numpy.array(MINFWHM_GSM),
     numpy.array(FACTOR), numpy.array(MINFWHM_GAUSSSLIT),
     numpy.array(FACTOR), numpy.array(MINFWHM_FULLGAUSS),
     numpy.array(FACTOR),numpy.array(MINFWHM_HYBRID),
     numpy.array(FACTOR)[-2:],numpy.array(MINFWHM_HYBRID)[-2:] * 0 + 9.13 * (18.8 / 65),
     numpy.array(FACTOR)[0:2],( numpy.array(FACTOR)[0:2] * beam_dimension_at_slit_in_um) * (28.8 / (65-35) ),
     xtitle="n", ytitle="WAIST FWHM [um]", ylog=0,
     marker=["o","o", "o", "o", "o", "None","None"],
     legend=["Wofry","WofryGSM", "WofryGaussianSlit", "WofryAllGauss","Hybrid","GeomOpt at source","GeomOpt at slit"])