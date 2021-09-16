import numpy
import xraylib




def transfocator_guess_configuration(focal_f_target,deltas=[0.999998],radii=[500e-4], verbose=1):

    nn = len(radii)

    ncombinations = 2**nn
    str1 = numpy.binary_repr(ncombinations-1)


    Farray = numpy.zeros(ncombinations)

    for i in range(ncombinations):
        fmt = "0%dd" % nn
        str1 = ("%" + fmt) % int(numpy.binary_repr(i))
        # print(">>> >>>", i , nn, str1)
        invF = 0
        for j in range(nn):
            if float(str1[j]) != 0:
                invF += 2 * deltas[j] / radii[j] * float(str1[j])
        if invF != 0:
            Farray[i] = 1.0 / invF
        else:
            Farray[i] = 1e5

    # print(Farray)


    iarg = numpy.argmin( numpy.abs(focal_f_target - Farray))

    if verbose:
        print(">>>> optimum: ", iarg )
        print(">>>> optimum for f=%g: " % focal_f_target, "%07d" % int(numpy.binary_repr(iarg)), Farray[iarg] )

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




"""

## TRANSFOCATOR @ 65
Transfocator 2D with 7 axis and 11 lenses
 - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
 - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
 - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
 - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
 - 1×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
 - 2×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
 - 4×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
## TRANSFOCATOR @ 65
Transfocator 1DH with 6 axis and 7 lenses
 - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
 - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
 - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
 - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
 - 1×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
 - 2×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)



## TRANSFOCATOR 2D @ 170
Transfocator 2D with 9 axis and 20 lenses
 - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
 - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
 - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
 - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
 - 1×Be lenses, r=300.0 μm, D=1.0 mm (2R_0=1079 μm)
 - 1×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
 - 2×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
 - 4×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
 - 8×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)

## TRANSFOCATOR 1DH @ 170
Transfocator with 4 axis and 4 lenses
 - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
 - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
 - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
 - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)


"""

if __name__ == "__main__":



    symbol = "Be"
    density = 1.845
    photon_energy_ev = 7000.0

    delta = 1.0 - xraylib.Refractive_Index_Re(symbol,photon_energy_ev*1e-3,density)

    print("delta: %g" % delta)

    # f1 in 15-85
    focal_f_target = 30.0


#
# ----------------------  TF 1 ---------------------------------
#
# - 1×Be lenses, r=5000.0 μm, D=1.0 mm (2R_0=4405 μm)
# - 1×Be lenses, r=2000.0 μm, D=1.0 mm (2R_0=2786 μm)
# - 1×Be lenses, r=1000.0 μm, D=1.0 mm (2R_0=1970 μm)
# - 1×Be lenses, r=500.0 μm, D=1.0 mm (2R_0=1393 μm)
# - 1×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
# - 2×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)
# - 4×Be lenses, r=200.0 μm, D=1.0 mm (2R_0=881 μm)


    fwanted = numpy.linspace(5,85,50)
    ffound = numpy.zeros_like(fwanted)

    for ii,focal_f_target in enumerate(fwanted):
        a = transfocator_guess_configuration(focal_f_target,
                                             deltas=[delta]*7,
                                             radii=[5000e-6,2000e-6,1000e-6,500e-6,200e-6,100e-5,50e-6])

        ffound[ii] = a

    print(ffound)
    # from srxraylib.plot.gol import plot, set_qt
    # set_qt()
    # plot(fwanted,ffound, xtitle="f wanted [m]", ytitle="f found [m]", title="TF1")


#
# ----------------------  TF 2 ---------------------------------
#
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


    ffound2 = numpy.zeros_like(fwanted)

    for ii,focal_f_target in enumerate(fwanted):
        a = transfocator_guess_configuration(focal_f_target,
                                             deltas=[delta]*9,
                                             radii=[5000e-6,2000e-6,1000e-6,500e-6,300e-6,200e-5,100e-6,50e-6,25e-6])

        ffound2[ii] = a

    print(ffound2)
    from srxraylib.plot.gol import plot, set_qt
    set_qt()
    plot(fwanted, fwanted,
         fwanted, ffound,
         fwanted,ffound2, xtitle="f wanted [m]", ytitle="f found [m]",
         legend=["ideal","TF1","TF2"],linestyle=[":",None,None])