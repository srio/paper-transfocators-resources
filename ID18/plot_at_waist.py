import numpy
from srxraylib.plot.gol import plot


use_real_lens = False

UP_TO_MODE        = [0,0,50,50]
USE_GAUSSIAN_SLIT = [True,False,True,False]

TMP_X = []
TMP_Y1 = []
TMP_Y2 = []
TMP_Y3 = []
TMP_Y4 = []
TMP_Y5 = []

for ii in range(len(UP_TO_MODE)):


        up_to_mode = UP_TO_MODE[ii]
        use_gaussian_slit = USE_GAUSSIAN_SLIT[ii]

        SIGMAS = ["0.1","0.2","0.5","1.0","1.5","2.0","4.0","6.0"]
        SIGMASF = []
        DISTANCE = []
        FWHM = []
        ICENTER = []
        ITOTAL = []
        LEGEND = []

        import scipy.constants as codata
        wavelength = codata.h*codata.c/codata.e/10000
        a = 234e-6 / 2
        p = 65.0
        print("N = ",a**2 / (wavelength * p))

        if use_gaussian_slit:
                gauss_add = "G"
        else:
                gauss_add = ""

        if use_real_lens:
                real_lens_add = "R"
        else:
                real_lens_add = ""

        if up_to_mode > 0:
                up_to_mode_add = "M"
        else:
                up_to_mode_add = ""

        for i in range(len(SIGMAS)):
                filein = "data_evolution/tmp%s%s%s%s.dat" % (up_to_mode_add,gauss_add, real_lens_add, SIGMAS[i])
                print(">>>>> ", filein)
                a1 = numpy.loadtxt(filein)
                print(a1.shape)
                distance1 = a1[:,0]
                fwhm1 = a1[:,1]
                itotal1 = a1[:,2]
                icenter1 = a1[:,3]
                DISTANCE.append(distance1)
                FWHM.append(fwhm1)
                ICENTER.append(icenter1)
                ITOTAL.append(itotal1)
                SIGMASF.append(float(SIGMAS[i]))
                slit_size_in_um = 125.0 / 2.35 * float(SIGMAS[i])
                s = slit_size_in_um * 1e-6 / 2
                pp = 35

                pa = 65 - pp
                # N2 =  p * s ** 2 / (wavelength * pp ** 2)
                N2 = p * s ** 2 / (wavelength * pa ** 2)
                print("Effect for slit aperture less than [um] = ", 2e6 * numpy.sqrt(pp ** 2 * wavelength / p))

                LEGEND.append(r'$a$=%s $\sigma_a$; N=%3.2f' % (SIGMAS[i], N2))



        WAISTPOSITION = []
        ICENTERATWAIST = []
        ITOTALATWAIST = []
        FWHMATWAIST = []
        FWHMAT99 = []

        for i in range(len(DISTANCE)):
                iMin1 = numpy.argmax(ICENTER[i])
                WAISTPOSITION.append(DISTANCE[i][iMin1])
                ITOTALATWAIST.append(ITOTAL[i][iMin1])
                ICENTERATWAIST.append(ICENTER[i][iMin1])
                FWHMATWAIST.append(FWHM[i][iMin1])
                FWHMAT99.append(
                        numpy.interp(
                        99.0,
                        numpy.array(DISTANCE[i]),
                        numpy.array(FWHM[i])))
                print("Minima found for: %g" % (WAISTPOSITION[i]))


        print("SIGMASF = ", SIGMASF)
        print("WAISTPOSITION = ", WAISTPOSITION)
        print("ICENTERATWAIST = ", ICENTERATWAIST)
        print("FWHMAT99 = ", FWHMAT99)
#
#         print("DISTANCE = ", DISTANCE[0])
#
#
        TMP_X.append(SIGMASF)
        TMP_Y1.append((WAISTPOSITION))
        TMP_Y2.append(ICENTERATWAIST)
        TMP_Y3.append(ITOTALATWAIST)
        TMP_Y4.append(FWHMATWAIST)
        TMP_Y5.append(FWHMAT99)


plot(   numpy.array(TMP_X[0]), numpy.array(TMP_Y1[0])/28.2,
        numpy.array(TMP_X[0]), numpy.array(TMP_Y1[1])/28.2,
        numpy.array(TMP_X[0]), numpy.array(TMP_Y1[2])/28.2,
        numpy.array(TMP_X[0]), numpy.array(TMP_Y1[3])/28.2,
        xlog=1, ylog=1, xtitle="n=a/(125/2.35)",ytitle="waist position over f",
        legend=["Gaussian slit","Rectangular slit","Gaussian slit Multimode","Rectangular slit Multimode",],
        linestyle=["--",None,'--',None],
        color=['red','red','blue','blue'],
        show=0)

plot(   numpy.array(TMP_X[0]), numpy.array(TMP_Y2[0]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y2[1]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y2[2]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y2[3]),
        xlog=1, ylog=1, xtitle="n=a/(125/2.35)",ytitle="Intensity at waist position",
        legend=["Gaussian slit","Rectangular slit","Gaussian slit Multimode","Rectangular slit Multimode",],
        linestyle=["--",None,'--',None],
        color=['red','red','blue','blue'],
        show=0)

plot(   numpy.array(TMP_X[0]), numpy.array(TMP_Y3[0]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y3[1]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y3[2]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y3[3]),
        xlog=1, ylog=0, xtitle="n=a/(125/2.35)",ytitle="Integrated Intensity at waist position",
        legend=["Gaussian slit","Rectangular slit","Gaussian slit Multimode","Rectangular slit Multimode",],
        linestyle=["--",None,'--',None],
        color=['red','red','blue','blue'],
        show=0)

plot(   numpy.array(TMP_X[0]), numpy.array(TMP_Y4[0]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y4[1]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y4[2]),
        numpy.array(TMP_X[0]), numpy.array(TMP_Y4[3]),
        xlog=0, ylog=0, xtitle="n=a/(125/2.35)",ytitle="FWHM at waist position [um]",
        legend=["Gaussian slit","Rectangular slit","Gaussian slit Multimode","Rectangular slit Multimode",],
        linestyle=["--",None,'--',None],
        color=['red','red','blue','blue'],
        show=0)

plot(   (125/2.35) * numpy.array(TMP_X[0]), numpy.array(TMP_Y5[0]),
        (125/2.35) * numpy.array(TMP_X[0]), numpy.array(TMP_Y5[1]),
        (125/2.35) * numpy.array(TMP_X[0]), numpy.array(TMP_Y5[2]),
        (125/2.35) * numpy.array(TMP_X[0]), numpy.array(TMP_Y5[3]),
        xlog=0, ylog=0, xtitle="a [um]",ytitle="FWHM 99 m from lens [um]",
        legend=["Gaussian slit","Rectangular slit","Gaussian slit Multimode","Rectangular slit Multimode",],
        linestyle=["--",None,'--',None],
        color=['red','red','blue','blue'],
        show=1)

