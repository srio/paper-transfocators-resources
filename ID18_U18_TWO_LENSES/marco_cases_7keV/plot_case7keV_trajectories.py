import numpy
from srxraylib.plot.gol import plot
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

    FILEROOT = ["trajectories_h/case7keV_2h","trajectories_v/case7keV_2v"]
    # fileroot = "trajectories_h/case7keV_2h"
    # fileroot = "trajectories_v/case7keV_2v"
    nfiles = 200

    #
    #
    #
    Index = numpy.arange(nfiles)

    for i,fileroot in enumerate(FILEROOT):

        if i == 0: # H
            marco_f1 = [41.69,29.73,100.00]
            marco_f2 = [26.15,24.27,25.16]
            marco_fwhm = [20.0, 30.0, 11.5]
            mc = numpy.loadtxt("../../GSM_MC/tmp_h.dat")
            slit = 40.0
        else:
            marco_f1 = [49.73,45.71,100]
            marco_f2 = [39.72,57.82,27.16]
            marco_fwhm = [20.0, 30.1, 7.5]
            mc = numpy.loadtxt("../../GSM_MC/tmp_v.dat")
            slit = 206.8

        F1 = []
        F2 = []
        F2theory1 = []
        F2theory2 = []
        R1 = []
        R2 = []
        FWHM = []
        M = []
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

            ff_source_at_id, mm_source_at_id = get_f2(f1=F1[index],
                                     position_source=0.0,  # source at source
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=True)

            ff_source_at_slit, mm_source_at_slit = get_f2(f1=F1[index],
                                     position_source=36.0, #0.0,  # source at slit
                                     position_lens1=65.0,
                                     position_lens2=170.0,
                                     position_sample=200.0,
                                     verbose=True)

            M.append(mm_source_at_slit)

            F2theory1.append(ff_source_at_id)
            F2theory2.append(ff_source_at_slit)


        print("F:", len(F1), len(F2))
        print("R:", len(R1), len(R2))

        from scipy.signal import savgol_filter

        F2theory1smooth = numpy.array(F2theory1)
        F2theory2smooth = numpy.array(F2theory2)
        Msmooth = numpy.array(M)

        for i in range(15):
            F2theory1smooth = savgol_filter(F2theory1smooth, 11, 1)  # window size 51, polynomial order 3
            F2theory2smooth = savgol_filter(F2theory2smooth, 11, 1)  # window size 51, polynomial order 3
            # Msmooth = savgol_filter(Msmooth, 11, 1)                  # window size 51, polynomial order 3


        plot(numpy.array(F1), numpy.array(F2),
            numpy.array(F1), F2theory1smooth,
            numpy.array(F1), F2theory2smooth,
            marco_f1, marco_f2,
            mc[:,0], mc[:,1],
            legend=["Wofry1D","Geometrical Optics smoothed (source at ID)","Geometrical Optics (source at slit)", "Marco Table", "Marco GSM"],
            marker=[None,None,None,'x',None],
            linestyle=['-','-','-',"None",'-'],
            yrange=[10,60],
            xtitle="F1 [m]", ytitle="F2 [m]", title="trajectories F %s" %fileroot)
        # plot(1e6 * numpy.array(R1), 1e6 * numpy.array(R2), xtitle="R1 [um]", ytitle="R2 [um]", title="trajectories R %s" %fileroot)
        plot(numpy.array(F1), numpy.array(FWHM),
             numpy.array(marco_f1), numpy.array(marco_fwhm),
             mc[:, 0], mc[:, 3],
             numpy.array(F1), Msmooth * slit,
             marker=[None,'x',None,None],
             linestyle=['-',"None",'-','-'],
             yrange=[0,numpy.array(FWHM).max()*1.1],
             legend=["Wofry1D","Marco Table","Marco GSM","Geometrical optics"],
             xtitle="F1 [m]", ytitle="FWHM [um]", title="Sizes %s" %fileroot)

