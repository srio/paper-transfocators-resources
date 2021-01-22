import numpy
import xraylib

"""
transfocator_id18 : transfocator for id18:
        It can:

            1) guess the lens configuration (number of lenses for each type) for a given photon energy
            and target image size. Use transfocator_compute_configuration() for this task

            2) for a given transfocator configuration, compute the main optical parameters
                (image size, focal distance, focal position and divergence).
                Use transfocator_compute_parameters() for this task

            3) Performs full ray tracing. Use id18_ray_tracing() for this task

        Note that for the optimization and parameters calculations the transfocator configuration is
        given in keywords. For ray tracing calculations many parameters of the transfocator are hard coded
        with the values of id18

        See main program for examples.

        Dependencies:
            Numpy
            xraylib (to compute refracion indices)
            Shadow (for ray tracing only)
            matplotlib (for some plots of ray=tracing)

        Side effects:
            When running ray tracing some files are created.

        MODIFICATION HISTORY:
           2015-03-25 srio@esrf.eu, written
           2020-11-06 srio@esrf.eu adapted to id18

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright__ = "ESRF, 2015"


def transfocator_compute_configuration(photon_energy_ev,
            s_target,
            symbol =["Be" ,"Be", "Be"],
            density=[1.845,1.845,1.845],
            nlenses_max   = [15,    3,      1],
            nlenses_radii = [500e-6,1000e-6,1500e-6], lens_diameter=500e-6,
            sigmaz=6.6e-6,
            alpha = 0.55,
            tf_p=164.0, tf_q=36.0, verbose=1 ):
    """
    Computes the optimum transfocator configuration for a given photon energy and target image size.

    All length units are m

    :param photon_energy_ev: the photon energy in eV
    :param s_target:       the target image size in m.
    :param symbol:         the chemical symbol of the lens material of each type. Default symbol=["Be","Be","Be"]
    :param density:        the density of each type of lens. Default: density=[1.845,1.845,1.845]
    :param nlenses_max:    the maximum allowed number of lenases for each type of lens. nlenses_max = [15,3,1]
    :param nlenses_radii:  the radii in m of each type of lens.
    :param lens_diameter:    the physical diameter (acceptance) in m of the lenses. If different for each type of lens,
                            consider the smaller one. Default: lens_diameter=0.0005
    :param sigmaz:         the sigma (standard deviation) of the source in m
    :param alpha:          an adjustable parameter in [0,1](see doc). Default: 0.55 (it is 0.76 for pure Gaussian beams)
    :param tf_p:           the distance source-transfocator in m
    :param tf_q:           the distance transfocator-image in m
    :param:verbose:        set to 1 for verbose text output
    :return:               a list with the number of lenses of each type.

    """
    if s_target < 2.35*sigmaz*tf_q/tf_p:
        print("Source size FWHM is: %f um"%(1e6*2.35*sigmaz))
        print("Maximum Demagnifications is: %f "%(tf_p/tf_q))
        print("Minimum possible size is: %f um"%(1e6*2.35*sigmaz*tf_q/tf_p))
        print("Error: redefine size")
        return None

    deltas = [(1.0 - xraylib.Refractive_Index_Re(symbol[i],photon_energy_ev*1e-3,density[i])) \
              for i in range(len(symbol))]

    focal_q_target = _tansfocator_guess_focal_position( s_target, p=tf_p, q=tf_q, sigmaz=sigmaz, alpha=alpha, \
                                           lens_diameter=lens_diameter,method=2)

    focal_f_target = 1.0 / (1.0/focal_q_target + 1.0/tf_p)
    div_q_target = alpha  * lens_diameter / focal_q_target

    #corrections for extreme cases
    source_demagnified = 2.35*sigmaz*focal_q_target/tf_p
    if source_demagnified > lens_diameter: source_demagnified = lens_diameter

    s_target_calc = numpy.sqrt( (div_q_target*(tf_q-focal_q_target))**2 + source_demagnified**2)

    nlenses_target = _transfocator_guess_configuration(focal_f_target,deltas=deltas,\
                                                   nlenses_max=nlenses_max,radii=nlenses_radii, )
    if verbose:
        print("transfocator_compute_configuration: focal_f_target: %f"%(focal_f_target))
        print("transfocator_compute_configuration: focal_q_target: %f m"%(focal_q_target))
        print("transfocator_compute_configuration: s_target: %f um"%(s_target_calc*1e6))
        print("transfocator_compute_configuration: nlenses_target: ",nlenses_target)

    return nlenses_target

def transfocator_compute_parameters(photon_energy_ev, nlenses_target,\
            symbol=["Be","Be","Be"], density=[1.845,1.845,1.845],\
            nlenses_max = [15,3,1], nlenses_radii = [500e-6,1000e-6,1500e-6], lens_diameter=0.0005, \
            sigmaz=6.6e-6, alpha = 0.55, \
            tf_p=164.0, tf_q=36.0 ):
    """
    Computes the parameters of the optical performances of a given transgocator configuration.

    returns a l

    All length units are m

    :param photon_energy_ev:
    :param nlenses_target: a list with the lens configuration, i.e. the number of lenses of each type.
    :param symbol:         the chemical symbol of the lens material of each type. Default symbol=["Be","Be","Be"]
    :param density:        the density of each type of lens. Default: density=[1.845,1.845,1.845]
    :param nlenses_radii:  the radii in m of each type of lens. Default: nlenses_radii = [500e-6,1000e-6,1500e-6]
    :param lens_diameter:    the physical diameter (acceptance) in m of the lenses. If different for each type of lens,
                            consider the smaller one. Default: lens_diameter=0.0005
    :param sigmaz:         the sigma (standard deviation) of the source in m
    :param alpha:          an adjustable parameter in [0,1](see doc). Default: 0.55 (it is 0.76 for pure Gaussian beams)
    :param tf_p:           the distance source-transfocator in m
    :param tf_q:           the distance transfocator-image in m
    :return:               a list with parameters (image_size, lens_focal_distance,
                           focal_position from transfocator center, divergence of beam after the transfocator)
    """

    deltas = [(1.0 - xraylib.Refractive_Index_Re(symbol[i],photon_energy_ev*1e-3,density[i])) \
              for i in range(len(symbol))]
    focal_f = _transfocator_calculate_focal_distance( deltas=deltas,\
                                                       nlenses=nlenses_target,radii=nlenses_radii)
    focal_q = 1.0 / (1.0/focal_f - 1.0/tf_p)
    div_q = alpha  * lens_diameter / focal_q
    #corrections
    source_demagnified = 2.35*sigmaz*focal_q/tf_p
    if source_demagnified > lens_diameter: source_demagnified = lens_diameter
    s_target = numpy.sqrt( (div_q*(tf_q-focal_q))**2 + (source_demagnified)**2 )
    return (s_target,focal_f,focal_q,div_q)


def transfocator_nlenses_to_slots(nlenses,nlenses_max=None):
    """
    converts the transfocator configuration from a list of the number of lenses of each type,
    into a list of active (1) or inactive (0) actuators for the slots.

    :param nlenses: the list with number of lenses (e.g., [5,2,0]
    :param nlenses_max: the maximum number of lenses of each type, usually powers of two minus one.
                        E.g. [15,3,1]
    :return: a list of on (1) and off (0) slots, e.g., [1, 0, 1, 0, 0, 1, 0]
            (first type: 1*1+0*2+1*4+0*8=5, second type: 0*1+1*2=2, third type: 0*1=0)
    """
    if nlenses_max == None:
        nlenses_max = nlenses

    ss = []
    for i,iopt in enumerate(nlenses):
        if iopt > nlenses_max[i]:
            print("Error: i:%d, nlenses: %d, nlenses_max: %d"%(i,iopt,nlenses_max[i]))
        ncharacters = len("{0:b}".format(nlenses_max[i]))

        si = list( ("{0:0%db}"%(ncharacters)).format(int(iopt)) )
        si.reverse()
        ss += si

    on_off = [int(i) for i in ss]
    #print("transfocator_nlenses_to_slots: nlenses_max: ",nlenses_max," nlenses: ",nlenses," slots: ",on_off)
    return on_off


def _transfocator_calculate_focal_distance(deltas=[0.999998],nlenses=[1],radii=[500e-6]):

    inverse_focal_distance = 0.0
    for i,nlensesi in enumerate(nlenses):
        if nlensesi > 0:
            focal_distance_i = radii[i] / (2.*nlensesi*deltas[i])
            inverse_focal_distance += 1.0/focal_distance_i
    if inverse_focal_distance == 0:
        return 99999999999999999999999999.
    else:
        return 1.0/inverse_focal_distance


def _tansfocator_guess_focal_position( s_target, p=164.0, q=36.0, sigmaz=6.6e-6, \
                                      alpha=0.66, lens_diameter=0.0005, method=2):
    x = 1e15
    if method == 1: # simple sum
        AA = 2.35*sigmaz/p
        BB = -(s_target + alpha * lens_diameter)
        CC = alpha*lens_diameter*q

        cc = numpy.roots([AA,BB,CC])
        x = cc[1]
        return x

    if method == 2: # sum in quadrature
        AA = ( (2.35*sigmaz)**2)/(p**2)
        BB = 0.0
        CC = alpha**2 * lens_diameter**2 - s_target**2
        DD = - 2.0 * alpha**2 * lens_diameter**2 * q
        EE = alpha**2 * lens_diameter**2 * q**2
        cc = numpy.roots([AA,BB,CC,DD,EE])
        for i,cci in enumerate(cc):
            if numpy.imag(cci) == 0:
                return numpy.real(cci)

    return x

def _transfocator_guess_configuration(focal_f_target,deltas=[0.999998],nlenses_max=[15],radii=[500e-6]):

    nn = len(nlenses_max)

    ncombinations = (1+nlenses_max[0]) * (1+nlenses_max[1]) * (1+nlenses_max[2])

    icombinations = 0
    aa = numpy.zeros((3,ncombinations),dtype=int)
    bb = numpy.zeros(ncombinations)
    for i0 in range(1+nlenses_max[0]):
        for i1 in range(1+nlenses_max[1]):
            for i2 in range(1+nlenses_max[2]):
                aa[0,icombinations] = i0
                aa[1,icombinations] = i1
                aa[2,icombinations] = i2
                bb[icombinations] = focal_f_target - _transfocator_calculate_focal_distance(deltas=deltas,nlenses=[i0,i1,i2],radii=radii)
                icombinations += 1
    bb1 = numpy.abs(bb)
    ibest = bb1.argmin()

    return (aa[:,ibest]).tolist()


def main():
    photon_energy_kev = 10.0   #  float(input("Enter photon energy in keV: "))
    s_target_um       = 10.0   #  float(input("Enter target focal dimension in microns: "))

    symbols       = ["Be",   "Be",    "Be"]
    densities     = [1.845,  1.845,   1.845]
    nlenses_max   = [15,     3,       1]
    nlenses_radii = [500e-6, 1000e-6, 1500e-6]
    lens_diameter = 0.0005

    alpha = 0.55

    sigmaz = 6.6e-6
    tf_p = 164.0
    tf_q = 36.00

    nlenses_optimum = transfocator_compute_configuration(
            photon_energy_kev * 1e3,
            s_target_um * 1e-6,
            symbol=symbols,
            density=densities,
            nlenses_max = nlenses_max,
            nlenses_radii = nlenses_radii,
            lens_diameter = lens_diameter,
            sigmaz = sigmaz,
            alpha = alpha,
            tf_p = tf_p,
            tf_q = tf_q,
            verbose=0 )

    print("Optimum lens configuration is: ",nlenses_optimum)

    if nlenses_optimum == None:
        return

    print("Activate slots: ",transfocator_nlenses_to_slots(nlenses_optimum,nlenses_max=nlenses_max))

    # this calculates the parameters (image size, etc) for a given lens configuration

    (size, f, q_f, div) = transfocator_compute_parameters(
            photon_energy_kev*1e3,
            nlenses_optimum,
            symbol=symbols,
            density=densities,
            nlenses_max=nlenses_max,
            nlenses_radii=nlenses_radii,
            lens_diameter=lens_diameter,
            sigmaz=sigmaz,
            alpha=alpha,
            tf_p=tf_p,
            tf_q=tf_q )

    print("For given configuration ",nlenses_optimum," we get: ")
    print("  size: %f um, focal length: %f m, focal distance: %f m, divergence: %f urad: "%(1e6*size, f, q_f, 1e6*div))

if __name__ == "__main__":
    main()