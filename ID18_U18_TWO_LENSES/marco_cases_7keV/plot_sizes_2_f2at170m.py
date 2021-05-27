import numpy
from srxraylib.plot.gol import plot, set_qt
from silx.io.specfile import SpecFile

def get_f2(f1=28.2,
           position_source=0.0,
           position_lens1=65.0,
           position_lens2= 170.0,
           position_sample=200.0,
           verbose=True):

    p1 = position_lens1 - position_source
    q1 = 1 / (1 / f1 - 1 / p1)


    p2 = position_lens2 - (p1 + q1)
    q2 = position_sample - position_lens2

    f2 = 1.0 / (1 / p2 + 1 / q2)

    if verbose:
        D = position_lens2 - position_lens1
        print("D: %g, q1+p2: %g" % (D, q1+p2))
        print("p1: %g" % p1)
        print("q1: %g" % q1)
        print("p2: %g" % p2)
        print("q2: %g" % q2)
        print("D: %g, Sum: %g" % (D, p1+q1+p2+q2))

    M = (q1 / q1) * (q2 / p2)
    return f2, M


if __name__ == "__main__":

    import xraylib

    FILEROOT = ["sizes_2h/case7keV_2h", "sizes_2v/case7keV_2v"]

    nfiles = 200

    #
    #
    #


    for i,fileroot in enumerate(FILEROOT):

        Index = numpy.arange(nfiles)
        if i == 0: # H
            marco_f1 = [41.69,29.73,100.00]
            marco_f2 = [26.15,24.27,25.16]
            marco_fwhm = [20.0, 30.0, 11.5]
            mc = numpy.loadtxt("../../GSM_MC/f1_f2_scans_170_h.dat")
            slit = 40.0
            sourcesize = 70.57
            prefix = "h"
        else:
            marco_f1 = [49.73,45.71,100]
            marco_f2 = [39.72,57.82,27.16]
            marco_fwhm = [20.0, 30.1, 7.5]
            mc = numpy.loadtxt("../../GSM_MC/f1_f2_scans_170_v.dat")
            slit = 206.8
            sourcesize = 15.02
            prefix = "v"




        a = numpy.loadtxt("trajectories_precalculated/f1_vs_f2_case2%s.dat" % prefix)
        F1 = a[:, 0].copy()
        F2 = a[:, 1].copy()




        #
        # read files with sizes
        #
        FWHM = []
        for index in Index:
            print(">>>>>>>>> opening file: ","%s_spectral_density_%03d.dat" % (fileroot, index))
            sf = SpecFile("%s_spectral_density_%03d.dat" % (fileroot, index))
            s1 = sf[0]
            fwhm = s1.scan_header_dict["UFWHM"]
            FWHM.append(float(fwhm))



        F2theory1 = []
        F2theory2 = []
        Msource_at_id = []
        Msource_at_slit = []

        for index in Index:
            ff_source_at_id, mm_source_at_id = get_f2(f1=F1[index],
                                     position_source=0.0,  # source at source
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=False)

            ff_source_at_slit, mm_source_at_slit = get_f2(f1=F1[index],
                                     position_source=36.0, #0.0,  # source at slit
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=False)

            Msource_at_id.append(mm_source_at_id)
            Msource_at_slit.append(mm_source_at_slit)

            F2theory1.append(ff_source_at_id)
            F2theory2.append(ff_source_at_slit)


        print("F:", len(F1), len(F2))

        from scipy.signal import savgol_filter

        F2theory1smooth = numpy.array(F2theory1)
        F2theory2smooth = numpy.array(F2theory2)
        Msource_at_id = numpy.array(Msource_at_id)
        Msource_at_slit = numpy.array(Msource_at_slit)

        # for i in range(15):
        #     F2theory1smooth = savgol_filter(F2theory1smooth, 11, 1)  # window size 51, polynomial order 3
        #     F2theory2smooth = savgol_filter(F2theory2smooth, 11, 1)  # window size 51, polynomial order 3
        #     # Msmooth = savgol_filter(Msmooth, 11, 1)                  # window size 51, polynomial order 3





        plot(numpy.array(F1), numpy.array(F2),
            numpy.array(F1), F2theory1smooth,
            numpy.array(F1), F2theory2smooth,
            marco_f1, marco_f2,
            mc[:,0], mc[:,1],
            legend=["Wofry1D","Geometrical Optics smoothed (source at ID)","Geometrical Optics (source at slit)", "Marco Table", "Marco GSM"],
            marker=[None,None,None,'x',None],
            linestyle=['-',':',':',"None",'-'],
            yrange=[10,60],
            xtitle="F1 [m]", ytitle="F2 [m]", title="trajectories F %s" %fileroot)
         plot(numpy.array(F1), numpy.array(FWHM),
             numpy.array(F1), Msource_at_id * sourcesize,
             numpy.array(F1), Msource_at_slit * slit,
             numpy.array(marco_f1), numpy.array(marco_fwhm),
             mc[:, 0], mc[:, 3],
             marker=[None,'.','.','x',None],
             linestyle=['-',"None","None","None",'-'],
             yrange=[0,numpy.array(FWHM).max()*1.1],
             legend=["Wofry1D","Geometrical optics (source at ID)","Geometrical optics (source at slit)","Marco Table","Marco GSM",],
             xtitle="F1 [m]", ytitle="FWHM [um]", title="Sizes %s" %fileroot)

