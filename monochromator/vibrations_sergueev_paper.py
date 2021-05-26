#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
from srxraylib.plot.gol import plot, plot_image


def run_source():
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()


    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 100000
    oe0.PH1 = 7000.0
    oe0.SIGDIX = 5e-16
    oe0.SIGDIZ = 5e-06
    oe0.SIGMAX = 5e-16
    oe0.SIGMAZ = 5e-06
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0



    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    return beam

def run_beamline(beam, x_rot1=2e-10, x_rot2=2e-10):
    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    # beam = Shadow.Beam()
    # oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()

    # #
    # # Define variables. See meaning of variables in:
    # #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    # #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    # #
    #
    # oe0.FDISTR = 3
    # oe0.F_PHOT = 0
    # oe0.HDIV1 = 0.0
    # oe0.HDIV2 = 0.0
    # oe0.IDO_VX = 0
    # oe0.IDO_VZ = 0
    # oe0.IDO_X_S = 0
    # oe0.IDO_Y_S = 0
    # oe0.IDO_Z_S = 0
    # oe0.ISTAR1 = 5676561
    # oe0.NPOINT = 100000
    # oe0.PH1 = 7000.0
    # oe0.SIGDIX = 5e-16
    # oe0.SIGDIZ = 5e-06
    # oe0.SIGMAX = 5e-16
    # oe0.SIGMAZ = 5e-06
    # oe0.VDIV1 = 0.0
    # oe0.VDIV2 = 0.0

    oe1.DUMMY = 100.0
    oe1.FILE_REFL = b'/users/srio/Oasys/Si5_35.111'
    oe1.FWRITE = 1
    oe1.F_CENTRAL = 1
    oe1.F_CRYSTAL = 1
    oe1.F_MOVE = 1
    oe1.PHOT_CENT = 7000.0
    oe1.R_LAMBDA = 5000.0
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 73.5910617052
    oe1.T_REFLECTION = 73.5910617052
    oe1.T_SOURCE = 47.0
    oe1.X_ROT = x_rot1 # = 2e-10

    oe2.ALPHA = 180.0
    oe2.DUMMY = 100.0
    oe2.FILE_REFL = b'/users/srio/Oasys/Si5_35.111'
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.F_MOVE = 1
    oe2.PHOT_CENT = 7000.0
    oe2.R_LAMBDA = 5000.0
    oe2.T_IMAGE = 3.0
    oe2.T_INCIDENCE = 73.5910617052
    oe2.T_REFLECTION = 73.5910617052
    oe2.T_SOURCE = 0.0
    oe2.X_ROT = x_rot2 # 3e-10

    oe3.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0006194038669230887, 0.0])
    oe3.DUMMY = 100.0
    oe3.FHIT_C = 1
    oe3.FILE_R_IND_IMA = b'/users/srio/Oasys/Be.dat'
    oe3.FMIRR = 4
    oe3.FSHAPE = 2
    oe3.FWRITE = 3
    oe3.F_EXT = 1
    oe3.F_REFRAC = 1
    oe3.F_R_IND = 2
    oe3.PARAM = 0.00030970193346154434
    oe3.RLEN2 = 0.0005
    oe3.RWIDX2 = 0.0005
    oe3.T_IMAGE = 2.5e-05
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = 180.0
    oe3.T_SOURCE = 0.0

    oe4.CCC = numpy.array([1.0, 1.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 0.0006194038669230887, 0.0])
    oe4.DUMMY = 100.0
    oe4.FHIT_C = 1
    oe4.FILE_R_IND_OBJ = b'/users/srio/Oasys/Be.dat'
    oe4.FMIRR = 4
    oe4.FSHAPE = 2
    oe4.FWRITE = 3
    oe4.F_CONVEX = 1
    oe4.F_EXT = 1
    oe4.F_REFRAC = 1
    oe4.F_R_IND = 1
    oe4.PARAM = 0.00030970193346154434
    oe4.RLEN2 = 0.0005
    oe4.RWIDX2 = 0.0005
    oe4.T_IMAGE = 40.0
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 2.5e-05

    # # Run SHADOW to create the source
    #
    # if iwrite:
    #     oe0.write("start.00")
    #
    # beam.genSource(oe0)
    #
    # if iwrite:
    #     oe0.write("end.00")
    #     beam.write("begin.dat")

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

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    return beam


def run_loop(sigma1_urad = 0.0, sigma2_urad = 0.0):


    sigma1_deg = sigma1_urad * 1e-6 * 180 / numpy.pi
    sigma2_deg = sigma2_urad * 1e-6 * 180 / numpy.pi

    r1 = numpy.random.normal(loc=0.0, scale=sigma1_deg, size=nruns)
    r2 = numpy.random.normal(loc=0.0, scale=sigma2_deg, size=nruns) # -r1 #

    sizes = numpy.zeros_like(distances)
    Sizes = numpy.zeros((nruns, distances.size))
    Centers = numpy.zeros((nruns, distances.size))
    RMSs = numpy.zeros((nruns, distances.size))

    HISTO_X = []
    HISTO_Y = []

    for i in range(nruns):
        beam = run_beamline(beam_source.duplicate(), x_rot1=r1[i], x_rot2=r2[i])

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
    beam_source = run_source()


    nruns = 500
    distances = numpy.linspace(-20, 10, 31)


    sizes0, hx0, hy0 = run_loop(sigma1_urad=0.0, sigma2_urad=0.0)
    sizes1, hx1, hy1 = run_loop(sigma1_urad=0.5, sigma2_urad=0.5)
    sizes2, hx2, hy2 = run_loop(sigma1_urad=1, sigma2_urad=1)




    plot(distances, 1e6 * sizes0,
         distances, 1e6 * sizes1,
         distances, 1e6 * sizes2,
         legend=["0 urad RMS", "0.5 urad RMS", "1 urad RMS"],
         xtitle="Distance from focal point [m]", ytitle="Beam size RMS")

    plot(numpy.array(hx0[0]), numpy.array(hy0[0]),
         numpy.array(hx1[0]), numpy.array(hy1[0]),
         numpy.array(hx2[0]), numpy.array(hy2[0]),
         )

    plot(numpy.array(hx0[distances.size // 2]), numpy.array(hy0[distances.size // 2]),
         numpy.array(hx1[distances.size // 2]), numpy.array(hy1[distances.size // 2]),
         numpy.array(hx2[distances.size // 2]), numpy.array(hy2[distances.size // 2]),
         )

    plot(numpy.array(hx0[-1]), numpy.array(hy0[-1]),
         numpy.array(hx1[-1]), numpy.array(hy1[-1]),
         numpy.array(hx2[-1]), numpy.array(hy2[-1]),
         )


    nx = distances.size
    ny = len(hx0[0])

    caustic0 = numpy.zeros((nx, ny))
    for i in range(nx):
        caustic0[i,:] = hy0[i]

    caustic1 = numpy.zeros((nx, ny))
    for i in range(nx):
        caustic1[i,:] = hy1[i]

    caustic2 = numpy.zeros((nx, ny))
    for i in range(nx):
        caustic2[i,:] = hy2[i]

    plot_image(caustic0, distances, numpy.array(hx0[0]), aspect='auto', title="0")
    plot_image(caustic1, distances, numpy.array(hx1[0]), aspect='auto', title="0.5")
    plot_image(caustic2, distances, numpy.array(hx2[0]), aspect='auto', title="1.0")