#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
from srxraylib.plot.gol import plot, plot_image
import pandas as pd


def run_beamline(beam, x_rot1=2e-10, x_rot2=2e-10, DCM_pos=60):
    """This function runs the beamline starting just after
    the hybrid widget, it takes the DCM crystal rotation misalignment 
    and position"""
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    oe0 = Shadow.OE()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()
    oe5 = Shadow.OE()
    oe6 = Shadow.OE()
    oe7 = Shadow.OE()
    oe8 = Shadow.OE()
    oe9 = Shadow.OE()
    oe10 = Shadow.OE()
    oe11 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.DUMMY = 100.0
    oe0.FILE_REFL = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/monochromator/Si_bragg_5_50.dat'
    oe0.FWRITE = 1
    oe0.F_CENTRAL = 1
    oe0.F_CRYSTAL = 1
    oe0.F_MOVE = 1
    oe0.PHOT_CENT = 7000.0
    oe0.R_LAMBDA = 5000.0
    oe0.T_IMAGE = 0.0
    oe0.T_INCIDENCE = 73.5910617052
    oe0.T_REFLECTION = 73.5910617052
    if DCM_pos == 38:
        oe0.T_SOURCE = 2.0
    elif DCM_pos == 60:
        oe0.T_SOURCE = 24.0
    else:
        raise RuntimeError(f'ERROR: Not DCM position has identified: {DCM_pos}')
    oe0.X_ROT = x_rot1  # = 2e-10 # pitch rotation

    oe1.DUMMY = 100.0
    oe1.FILE_REFL = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/monochromator/Si_bragg_5_50.dat'
    oe1.FWRITE = 1
    oe1.F_CENTRAL = 1
    oe1.F_CRYSTAL = 1
    oe1.F_MOVE = 1
    oe1.PHOT_CENT = 7000.0
    oe1.R_LAMBDA = 5000.0
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 73.5910617052
    oe1.T_REFLECTION = 73.5910617052
    oe1.T_SOURCE = 0.0
    oe1.X_ROT = x_rot2 # pitch rotation

    oe2.DUMMY = 100.0
    oe2.FWRITE = 3
    oe2.F_REFRAC = 2
    oe2.F_SCREEN = 1
    oe2.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe2.K_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe2.N_SCREEN = 1
    oe2.RX_SLIT = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe2.RZ_SLIT = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 0.0
    oe2.T_REFLECTION = 180.0
    if DCM_pos == 38:
        oe2.T_SOURCE = 27.0 #For DCM @ 38 m
    elif DCM_pos == 60:
        oe2.T_SOURCE = 5.0 #For DCM @ 60 m
    else:
        raise RuntimeError(f'ERROR: Not DCM position has identified: {DCM_pos}')

    oe3.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0011606, 0.0])
    oe3.CIL_ANG = 90.0
    oe3.DUMMY = 100.0
    oe3.FCYL = 1
    oe3.FHIT_C = 1
    oe3.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe3.FMIRR = 4
    oe3.FSHAPE = 2
    oe3.FWRITE = 3
    oe3.F_EXT = 1
    oe3.F_REFRAC = 1
    oe3.F_R_IND = 2
    oe3.PARAM = 0.0005803
    oe3.RLEN2 = 0.0005
    oe3.RWIDX2 = 0.0005
    oe3.T_IMAGE = 5e-06
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = 180.0
    oe3.T_SOURCE = 0.0

    oe4.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 0.0011606, 0.0])
    oe4.CIL_ANG = 90.0
    oe4.DUMMY = 100.0
    oe4.FCYL = 1
    oe4.FHIT_C = 1
    oe4.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe4.FMIRR = 4
    oe4.FSHAPE = 2
    oe4.FWRITE = 3
    oe4.F_CONVEX = 1
    oe4.F_EXT = 1
    oe4.F_REFRAC = 1
    oe4.F_R_IND = 1
    oe4.PARAM = 0.0005803
    oe4.RLEN2 = 0.0005
    oe4.RWIDX2 = 0.0005
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 5e-06

    oe5.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0013842, 0.0])
    oe5.DUMMY = 100.0
    oe5.FCYL = 1
    oe5.FHIT_C = 1
    oe5.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe5.FMIRR = 4
    oe5.FSHAPE = 2
    oe5.FWRITE = 3
    oe5.F_EXT = 1
    oe5.F_REFRAC = 1
    oe5.F_R_IND = 2
    oe5.PARAM = 0.0006921
    oe5.RLEN2 = 0.0005
    oe5.RWIDX2 = 0.0005
    oe5.T_IMAGE = 5e-06
    oe5.T_INCIDENCE = 0.0
    oe5.T_REFLECTION = 180.0
    oe5.T_SOURCE = 0.0

    oe6.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 0.0013842, 0.0])
    oe6.DUMMY = 100.0
    oe6.FCYL = 1
    oe6.FHIT_C = 1
    oe6.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe6.FMIRR = 4
    oe6.FSHAPE = 2
    oe6.FWRITE = 3
    oe6.F_CONVEX = 1
    oe6.F_EXT = 1
    oe6.F_REFRAC = 1
    oe6.F_R_IND = 1
    oe6.PARAM = 0.0006921
    oe6.RLEN2 = 0.0005
    oe6.RWIDX2 = 0.0005
    oe6.T_IMAGE = 0.0
    oe6.T_INCIDENCE = 0.0
    oe6.T_REFLECTION = 180.0
    oe6.T_SOURCE = 5e-06

    oe7.DUMMY = 100.0
    oe7.FWRITE = 3
    oe7.F_REFRAC = 2
    oe7.F_SCREEN = 1
    oe7.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe7.K_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe7.N_SCREEN = 1
    oe7.RX_SLIT = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe7.RZ_SLIT = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe7.T_IMAGE = 0.0
    oe7.T_INCIDENCE = 0.0
    oe7.T_REFLECTION = 180.0
    oe7.T_SOURCE = 105.0

    oe8.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.000728, 0.0])
    oe8.CIL_ANG = 90.0
    oe8.DUMMY = 100.0
    oe8.FCYL = 1
    oe8.FHIT_C = 1
    oe8.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe8.FMIRR = 4
    oe8.FSHAPE = 2
    oe8.FWRITE = 3
    oe8.F_EXT = 1
    oe8.F_REFRAC = 1
    oe8.F_R_IND = 2
    oe8.PARAM = 0.000364
    oe8.RLEN2 = 0.0005
    oe8.RWIDX2 = 0.0005
    oe8.T_IMAGE = 5e-06
    oe8.T_INCIDENCE = 0.0
    oe8.T_REFLECTION = 180.0
    oe8.T_SOURCE = 0.0

    oe9.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 0.000728, 0.0])
    oe9.CIL_ANG = 90.0
    oe9.DUMMY = 100.0
    oe9.FCYL = 1
    oe9.FHIT_C = 1
    oe9.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe9.FMIRR = 4
    oe9.FSHAPE = 2
    oe9.FWRITE = 3
    oe9.F_CONVEX = 1
    oe9.F_EXT = 1
    oe9.F_REFRAC = 1
    oe9.F_R_IND = 1
    oe9.PARAM = 0.000364
    oe9.RLEN2 = 0.0005
    oe9.RWIDX2 = 0.0005
    oe9.T_IMAGE = 0.0
    oe9.T_INCIDENCE = 0.0
    oe9.T_REFLECTION = 180.0
    oe9.T_SOURCE = 5e-06

    oe10.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0011058, 0.0])
    oe10.DUMMY = 100.0
    oe10.FCYL = 1
    oe10.FHIT_C = 1
    oe10.FILE_R_IND_IMA = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe10.FMIRR = 4
    oe10.FSHAPE = 2
    oe10.FWRITE = 3
    oe10.F_EXT = 1
    oe10.F_REFRAC = 1
    oe10.F_R_IND = 2
    oe10.PARAM = 0.0005529
    oe10.RLEN2 = 0.0005
    oe10.RWIDX2 = 0.0005
    oe10.T_IMAGE = 5e-06
    oe10.T_INCIDENCE = 0.0
    oe10.T_REFLECTION = 180.0
    oe10.T_SOURCE = 0.0

    oe11.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 0.0011058, 0.0])
    oe11.DUMMY = 100.0
    oe11.FCYL = 1
    oe11.FHIT_C = 1
    oe11.FILE_R_IND_OBJ = b'/mntdirect/_users/reyesher/OASYS/EBSL1-ID18/Be_lens.dat'
    oe11.FMIRR = 4
    oe11.FSHAPE = 2
    oe11.FWRITE = 3
    oe11.F_CONVEX = 1
    oe11.F_EXT = 1
    oe11.F_REFRAC = 1
    oe11.F_R_IND = 1
    oe11.PARAM = 0.0005529
    oe11.RLEN2 = 0.0005
    oe11.RWIDX2 = 0.0005
    oe11.T_IMAGE = 30.0
    oe11.T_INCIDENCE = 0.0
    oe11.T_REFLECTION = 180.0
    oe11.T_SOURCE = 5e-06

    # run optical element 0

    print("    Running optical element: %d" % (0))
    if iwrite:
        oe0.write("start.00")

    beam.traceOE(oe0, 0)

    if iwrite:
        oe0.write("end.00")
        beam.write("star.00")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    #
    # run optical element 4
    #
    print("    Running optical element: %d" % (4))
    if iwrite:
        oe4.write("start.04")

    beam.traceOE(oe4, 4)

    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")

    #
    # run optical element 5
    #
    print("    Running optical element: %d" % (5))
    if iwrite:
        oe5.write("start.05")

    beam.traceOE(oe5, 5)

    if iwrite:
        oe5.write("end.05")
        beam.write("star.05")

    #
    # run optical element 6
    #
    print("    Running optical element: %d" % (6))
    if iwrite:
        oe6.write("start.06")

    beam.traceOE(oe6, 6)

    if iwrite:
        oe6.write("end.06")
        beam.write("star.06")

    #
    # run optical element 7
    #
    print("    Running optical element: %d" % (7))
    if iwrite:
        oe7.write("start.07")

    beam.traceOE(oe7, 7)

    if iwrite:
        oe7.write("end.07")
        beam.write("star.07")

    #
    # run optical element 8
    #
    print("    Running optical element: %d" % (8))
    if iwrite:
        oe8.write("start.08")

    beam.traceOE(oe8, 8)

    if iwrite:
        oe8.write("end.08")
        beam.write("star.08")

    #
    # run optical element 9
    #
    print("    Running optical element: %d" % (9))
    if iwrite:
        oe9.write("start.09")

    beam.traceOE(oe9, 9)

    if iwrite:
        oe9.write("end.09")
        beam.write("star.09")

    #
    # run optical element 10
    #
    print("    Running optical element: %d" % (10))
    if iwrite:
        oe10.write("start.10")

    beam.traceOE(oe10, 10)

    if iwrite:
        oe10.write("end.10")
        beam.write("star.10")

    #
    # run optical element 11
    #
    print("    Running optical element: %d" % (11))
    if iwrite:
        oe11.write("start.11")

    beam.traceOE(oe11, 11)

    if iwrite:
        oe11.write("end.11")
        beam.write("star.11")


    return beam


def run_loop(sigma1_urad = 0.0, sigma2_urad = 0.0, DCM_pos=38):


    sigma1_deg = sigma1_urad * 1e-6 * 180 / numpy.pi
    sigma2_deg = sigma2_urad * 1e-6 * 180 / numpy.pi
 
    
    #random normal distributions
    r1 = numpy.random.normal(loc=0.0, scale=sigma1_deg, size=nruns)
    r2 = numpy.random.normal(loc=0.0, scale=sigma2_deg, size=nruns) # -r1 #

    sizes = numpy.zeros_like(distances)
    Sizes = numpy.zeros((nruns, distances.size))
    Centers = numpy.zeros((nruns, distances.size))
    RMSs = numpy.zeros((nruns, distances.size))

    HISTO_X = []
    HISTO_Y = []

    for i in range(nruns):
        beam = run_beamline(beam_source.duplicate(), x_rot1=r1[i], x_rot2=r2[i], DCM_pos=DCM_pos)

        for j,distance in enumerate(distances):
            tmp = beam.duplicate()
            tmp.retrace(distance)

            tkt = tmp.histo1(3, xrange=[-300e-6,300e-6], nbins=150, nolost=1, ref=0, write=None, factor=1.0, calculate_widths=1,
                   calculate_hew=0)
            fwhm = tkt["fwhm"]
            if i == 0:
                HISTO_X.append(tkt["bin_center"])
                HISTO_Y.append(tkt["histogram"])
            else:
                HISTO_X[j] = tkt["bin_center"]
                HISTO_Y[j] += tkt["histogram"]

            # plot(tkt["bin_center"], tkt["histogram"], title="%d %d " % (i,j) )
            col3 = tmp.getshonecol(3, nolost=1)
            weight = numpy.ones_like(col3) #  = tmp.getshonecol(23, nolost=1)
            sdev = tmp.get_standard_deviation(3, nolost=1, ref=0)
            center = numpy.average(col3, weights=weight)

            vv = numpy.average((col3) ** 2, weights=weight)
            rms = numpy.sqrt(vv)

            Sizes[i,j] = sdev
            Centers[i,j] = center
            RMSs[i,j]  = rms


    return numpy.abs(RMSs).sum(axis=0) / nruns, HISTO_X, HISTO_Y


if __name__  == "__main__":

    beam_source = Shadow.Beam()
    beam_source.load('Hy_screen_beam_1M.dat')


    nruns = 300    
    distances = numpy.linspace(-15, 15, 31)
    #@38 m or 60 m from the source
    DCM_pos = 38 #60    

    sizes0, hx0, hy0 = run_loop(sigma1_urad=0.0, sigma2_urad=0.0, DCM_pos=DCM_pos)
    sizes1, hx1, hy1 = run_loop(sigma1_urad=0.5, sigma2_urad=0.5, DCM_pos=DCM_pos)    
    sizes2, hx2, hy2 = run_loop(sigma1_urad=1, sigma2_urad=1, DCM_pos=DCM_pos)
    
    # saving into a csv
    data = {'distances':distances,'0 urad RMS':sizes0,'0.5 urad RMS':sizes1,'1 urad RMS':sizes2}
    df = pd.DataFrame(data)
    df.to_csv(f'DCM_vibs_dist_to_focal_point_pos_{DCM_pos}.csv')

    # plotting options
    # plot(distances, 1e6 * sizes0,
         # distances, 1e6 * sizes1,
         # distances, 1e6 * sizes2,
         # legend=["0 urad RMS", "0.5 urad RMS", "1 urad RMS"],
         # xtitle="Distance from focal point [m]", ytitle="Beam size RMS")

    # plot(numpy.array(hx0[0]), numpy.array(hy0[0]),
         # numpy.array(hx1[0]), numpy.array(hy1[0]),
         # numpy.array(hx2[0]), numpy.array(hy2[0]),
         # )

    # plot(numpy.array(hx0[distances.size // 2]), numpy.array(hy0[distances.size // 2]),
         # numpy.array(hx1[distances.size // 2]), numpy.array(hy1[distances.size // 2]),
         # numpy.array(hx2[distances.size // 2]), numpy.array(hy2[distances.size // 2]),
         # )

    # plot(numpy.array(hx0[-1]), numpy.array(hy0[-1]),
         # numpy.array(hx1[-1]), numpy.array(hy1[-1]),
         # numpy.array(hx2[-1]), numpy.array(hy2[-1]),
         # )


    # nx = distances.size
    # ny = len(hx0[0])

    # caustic0 = numpy.zeros((nx, ny))
    # for i in range(nx):
        # caustic0[i,:] = hy0[i]

    # caustic1 = numpy.zeros((nx, ny))
    # for i in range(nx):
        # caustic1[i,:] = hy1[i]

    # caustic2 = numpy.zeros((nx, ny))
    # for i in range(nx):
        # caustic2[i,:] = hy2[i]

    # plot_image(caustic0, distances, numpy.array(hx0[0]), aspect='auto', title="0")
    # plot_image(caustic1, distances, numpy.array(hx1[0]), aspect='auto', title="0.5")
    # plot_image(caustic2, distances, numpy.array(hx2[0]), aspect='auto', title="1.0")
