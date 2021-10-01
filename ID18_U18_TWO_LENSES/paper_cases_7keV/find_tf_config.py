import numpy
import xraylib




def transfocator_guess_configuration(focal_f_target, deltas=[0.999998], radii=[500e-4],
                                     initial_focal_distance=None, verbose=0):

    nn = len(radii)
    ncombinations = 2**nn
    Farray = numpy.zeros(ncombinations)
    # Rarray = numpy.zeros(ncombinations)

    for i in range(ncombinations):
        str1 = numpy.binary_repr(i, width=nn)
        if initial_focal_distance is None:
            invF = 0
            # invR = 0
        else:
            invF = 1.0 / initial_focal_distance
            # invR = invF / (2 * deltas[0])

        for j in range(nn):
            # if float(str1[j]) != 0:
            invF += 2 * deltas[j] / radii[j] * float(str1[j])
            # invR += 1 / radii[j] * float(str1[j])
            # try:
            #     print(">>>", i, nn, j, str1, float(str1[j]), 1/invF, 1e6/invR, ">", 1e6*radii[j])
            # except:
            #     print(">>>>>", i, nn, j, str1, float(str1[j]))
        if invF != 0:
            Farray[i] = 1.0 / invF
        else:
            Farray[i] = 1e5
        # if invR != 0:
        #     Rarray[i] = 1.0 / invR
        # else:
        #     Rarray[i] = 1e15
    # print(Farray)


    iarg = numpy.argmin( numpy.abs(focal_f_target - Farray))

    if verbose:
        # print(">>>> optimum: ", iarg )
        print(">>>> optimum for f=%g (idx %d): " % (focal_f_target, iarg), numpy.binary_repr(iarg, width=nn), Farray[iarg] )
        print("     cumulated radius: wanted R=%g found R=%g: " % (
            1e6*_transfocator_calculate_radius(delta=deltas[0], focal_distance=focal_f_target),
            1e6*_transfocator_calculate_radius(delta=deltas[0], focal_distance=Farray[iarg]),
            # 1e6*Rarray[iarg]
                                           ))
        print("     Initial focal distance: ", initial_focal_distance)

    return Farray[iarg]

def _transfocator_calculate_focal_distance(deltas=[0.999998],nlenses=[1],radii=[500e-4]):

    inverse_focal_distance = 0.0
    for i,nlensesi in enumerate(nlenses):
        if nlensesi > 0:
            focal_distance_i = radii[i] / (2.*nlensesi*deltas[i])
            inverse_focal_distance += 1.0/focal_distance_i
    if inverse_focal_distance == 0:
        return 99999999999999999999999999.
    else:
        return 1.0/inverse_focal_distance

def _transfocator_calculate_radius(delta=0.999998,focal_distance=10):
    radius = focal_distance * (2.*delta)
    return radius

if __name__ == "__main__":



    symbol = "Be"
    density = 1.845
    photon_energy_ev = 7000.0
    delta = 1.0 - xraylib.Refractive_Index_Re(symbol,photon_energy_ev*1e-3,density)
    print("delta: %g" % delta)

    # f1 in 15-85
    # focal_f_target = 30.0

    fwanted = numpy.linspace(2,85,50)
    # fwanted = numpy.array([2])


    # fpaper_tf1v = numpy.array([15.0, 42.2, 85.2, 42.2])
    # fpaper_tf1h = numpy.array([46.1, 25.1, 46.1, 25.1])
    #
    # fpaper_tf2v = numpy.array([22.2, 55.6, 27.8, 55.7])
    # fpaper_tf2h = numpy.array([26.5, 21.3, 31.8, 20.7])

    fpaper_tf1v = numpy.array([ 42.2 ])
    fpaper_tf1h = numpy.array([ 25.1 ])
    fpaper_tf2v = numpy.array([ 55.7 ])
    fpaper_tf2h = numpy.array([ 20.7 ])


# ## TRANSFOCATOR @ 65
# Transfocator 2D with 7 axis and 11 lenses
#  - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
#  - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
#  - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
#  - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
#  - 1×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
#  - 2×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
#  - 4×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
# ## TRANSFOCATOR @ 65
# Transfocator 1DH with 6 axis and 7 lenses
#  - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
#  - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
#  - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
#  - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
#  - 1×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
#  - 2×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)

    radii_tf1v = [5000e-6, 2000e-6, 1000e-6, 500e-6, 200e-6, 200e-6/2, 200e-6/4]
    radii_tf1h = [5000e-6, 2000e-6, 1000e-6, 500e-6, 200e-6, 200e-6/2]


# ## TRANSFOCATOR 2D @ 170
# Transfocator 2D with 9 axis and 20 lenses
#  - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
#  - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
#  - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
#  - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
#  - 1×Be lenses, r=300.0 μm, D=1.0 mm (2R_0=1079 μm)
#  - 1×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
#  - 2×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
#  - 4×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
#  - 8×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
#
# ## TRANSFOCATOR 1DH @ 170
# Transfocator with 4 axis and 4 lenses
#  - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
#  - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
#  - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
#  - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
    radii_tf2v = [5000e-6, 2000e-6, 1000e-6, 500e-6, 300e-6, 200e-6, 200e-6/2, 200e-6/4, 200e-6/8]
    radii_tf2h = [5000e-6, 2000e-6, 1000e-6, 500e-6]


    if False:
        ffound = numpy.zeros_like(fwanted)
        for ii,focal_f_target in enumerate(fwanted):
            a = transfocator_guess_configuration(focal_f_target, deltas=[delta]*len(radii_tf1v), radii=radii_tf1v)
            ffound[ii] = a
        # print(ffound)

        ffound2 = numpy.zeros_like(fwanted)
        for ii,focal_f_target in enumerate(fwanted):
            a = transfocator_guess_configuration(focal_f_target,
                                                 deltas=[delta]*len(radii_tf2v), radii=radii_tf2v)
            ffound2[ii] = a
        # print(ffound2)



        #
        # plot
        #
        from srxraylib.plot.gol import plot, set_qt
        set_qt()
        plot(fwanted, fwanted,
             fwanted, ffound,
             fwanted,ffound2,
             fpaper_tf1v, fpaper_tf1v,
             fpaper_tf2v, fpaper_tf2v,
             fpaper_tf1h, fpaper_tf1h,
             fpaper_tf2h, fpaper_tf2h,
             xtitle="f wanted [m]", ytitle="f found [m]",
             legend=["ideal","TF1","TF2","f wanted TF1 V","f wanted TF2 V","f wanted TF1 H","f wanted TF2 H"],
             linestyle=[":",None,None,"","","",""],
             marker=[None,None,None,'+','+','x','x'],
             title="2D focusing")


    #
    # TF1
    #
    fwanted_2d = numpy.zeros_like(fpaper_tf1h)
    ffound_2d = numpy.zeros_like(fpaper_tf1h)
    for i in range(fwanted_2d.size):
        fwanted_2d[i] = numpy.max( (fpaper_tf1h[i], fpaper_tf1v[i]))
        tmp = transfocator_guess_configuration(fwanted_2d[i], deltas=[delta]*len(radii_tf1v),
                                               radii=radii_tf1v, verbose=1)
        ffound_2d[i] = tmp


    fwanted_1d = numpy.zeros_like(fpaper_tf1h)
    ffound_1d = numpy.zeros_like(fpaper_tf1h)
    for i in range(fwanted_1d.size):
        fwanted_1d[i] = fpaper_tf1h[i]
        tmp = transfocator_guess_configuration(fwanted_1d[i], deltas=[delta]*len(radii_tf1h),
                                               radii=radii_tf1h, verbose=1, initial_focal_distance=ffound_2d[i])
        ffound_1d[i] = tmp

    print("TF1 V 2D f wanted, f found: ", fpaper_tf1v,ffound_2d)
    print("TF1 H 1D f wanted, f found: ", fpaper_tf1h,ffound_1d)


    #
    # TF2
    #
    fwanted_2d = numpy.zeros_like(fpaper_tf2h)
    ffound_2d = numpy.zeros_like( fpaper_tf2h)
    for i in range(fwanted_2d.size):
        fwanted_2d[i] = numpy.max( (fpaper_tf2h[i], fpaper_tf2v[i]))
        tmp = transfocator_guess_configuration(fwanted_2d[i], deltas=[delta]*len(radii_tf2v),
                                               radii=radii_tf2v, verbose=1)
        ffound_2d[i] = tmp


    fwanted_1d = numpy.zeros_like(fpaper_tf2h)
    ffound_1d  = numpy.zeros_like(fpaper_tf2h)
    for i in range(fwanted_1d.size):
        fwanted_1d[i] = fpaper_tf2h[i]
        tmp = transfocator_guess_configuration(fwanted_1d[i], deltas=[delta]*len(radii_tf2h),
                                               radii=radii_tf2h, verbose=1, initial_focal_distance=ffound_2d[i])
        ffound_1d[i] = tmp

    print("TF2 V 2D f wanted, f found: ", fpaper_tf2v,ffound_2d)
    print("TF2 H 1D f wanted, f found: ", fpaper_tf2h,ffound_1d)

    print(1.0 / (1/5000 + 1/1000 + 1/5000 + 1/500))
