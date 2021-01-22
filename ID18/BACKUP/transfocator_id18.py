import numpy
import xraylib
from transfocator_id18_lighted import transfocator_compute_configuration, transfocator_nlenses_to_slots, \
    transfocator_compute_parameters


"""
transfocator_id18 : transfocator for id18:
        It can:

            1) guess the lens configuration (number of lenses for each type) for a given photon energy
            and target image size. Use transfocator_compute_configuration() for this task

            2) for a given transfocator configuration, compute the main optical parameters
                (image size, focal distance, focal position and divergence).
                Use transfocator_compute_parameters() for this task

            3) Performs full ray tracing. Use id30b_ray_tracing() for this task

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

#
# def transfocator_compute_configuration(photon_energy_ev,s_target,\
#             symbol=["Be","Be","Be"], density=[1.845,1.845,1.845],\
#             nlenses_max = [15,3,1], nlenses_radii = [500e-4,1000e-4,1500e-4], lens_diameter=0.05, \
#             sigmaz=6.6e-4, alpha = 0.55, \
#             tf_p=6500, tf_q=13500, verbose=1 ):
#     """
#     Computes the optimum transfocator configuration for a given photon energy and target image size.
#
#     All length units are cm
#
#     :param photon_energy_ev: the photon energy in eV
#     :param s_target:       the target image size in cm.
#     :param symbol:         the chemical symbol of the lens material of each type. Default symbol=["Be","Be","Be"]
#     :param density:        the density of each type of lens. Default: density=[1.845,1.845,1.845]
#     :param nlenses_max:    the maximum allowed number of lenases for each type of lens. nlenses_max = [15,3,1]
#     :param nlenses_radii:  the radii in cm of each type of lens. Default: nlenses_radii = [500e-4,1000e-4,1500e-4]
#     :param lens_diameter:    the physical diameter (acceptance) in cm of the lenses. If different for each type of lens,
#                             consider the smaller one. Default: lens_diameter=0.05
#     :param sigmaz:         the sigma (standard deviation) of the source in cm
#     :param alpha:          an adjustable parameter in [0,1](see doc). Default: 0.55 (it is 0.76 for pure Gaussian beams)
#     :param tf_p:           the distance source-transfocator in cm
#     :param tf_q:           the distance transfocator-image in cm
#     :param:verbose:        set to 1 for verbose text output
#     :return:               a list with the number of lenses of each type.
#
#     """
#     if s_target < 2.35*sigmaz*tf_q/tf_p:
#         print("Source size FWHM is: %f um"%(1e4*2.35*sigmaz))
#         print("Maximum Demagnifications is: %f um"%(tf_p/tf_q))
#         print("Minimum possible size is: %f um"%(1e4*2.35*sigmaz*tf_q/tf_p))
#         print("Error: redefine size")
#         return None
#
#     deltas = [(1.0 - xraylib.Refractive_Index_Re(symbol[i],photon_energy_ev*1e-3,density[i])) \
#               for i in range(len(symbol))]
#
#     focal_q_target = _tansfocator_guess_focal_position( s_target, p=tf_p, q=tf_q, sigmaz=sigmaz, alpha=alpha, \
#                                            lens_diameter=lens_diameter,method=2)
#
#     focal_f_target = 1.0 / (1.0/focal_q_target + 1.0/tf_p)
#     div_q_target = alpha  * lens_diameter / focal_q_target
#
#     #corrections for extreme cases
#     source_demagnified = 2.35*sigmaz*focal_q_target/tf_p
#     if source_demagnified > lens_diameter: source_demagnified = lens_diameter
#
#     s_target_calc = numpy.sqrt( (div_q_target*(tf_q-focal_q_target))**2 + source_demagnified**2)
#
#     nlenses_target = _transfocator_guess_configuration(focal_f_target,deltas=deltas,\
#                                                    nlenses_max=nlenses_max,radii=nlenses_radii, )
#     if verbose:
#         print("transfocator_compute_configuration: focal_f_target: %f"%(focal_f_target))
#         print("transfocator_compute_configuration: focal_q_target: %f cm"%(focal_q_target))
#         print("transfocator_compute_configuration: s_target: %f um"%(s_target_calc*1e4))
#         print("transfocator_compute_configuration: nlenses_target: ",nlenses_target)
#
#     return nlenses_target
#
# def transfocator_compute_parameters(photon_energy_ev, nlenses_target,\
#             symbol=["Be","Be","Be"], density=[1.845,1.845,1.845],\
#             nlenses_max = [15,3,1], nlenses_radii = [500e-4,1000e-4,1500e-4], lens_diameter=0.05, \
#                                     sigmaz=6.6e-4, alpha=0.55, \
#                                     tf_p=6500, tf_q=13500 ):
#     """
#     Computes the parameters of the optical performances of a given transgocator configuration.
#
#     returns a l
#
#     All length units are cm
#
#     :param photon_energy_ev:
#     :param nlenses_target: a list with the lens configuration, i.e. the number of lenses of each type.
#     :param symbol:         the chemical symbol of the lens material of each type. Default symbol=["Be","Be","Be"]
#     :param density:        the density of each type of lens. Default: density=[1.845,1.845,1.845]
#     :param nlenses_max:    the maximum allowed number of lenases for each type of lens. nlenses_max = [15,3,1]
#                            TODO: remove (not used)
#     :param nlenses_radii:  the radii in cm of each type of lens. Default: nlenses_radii = [500e-4,1000e-4,1500e-4]
#     :param lens_diameter:    the physical diameter (acceptance) in cm of the lenses. If different for each type of lens,
#                             consider the smaller one. Default: lens_diameter=0.05
#     :param sigmaz:         the sigma (standard deviation) of the source in cm
#     :param alpha:          an adjustable parameter in [0,1](see doc). Default: 0.55 (it is 0.76 for pure Gaussian beams)
#     :param tf_p:           the distance source-transfocator in cm
#     :param tf_q:           the distance transfocator-image in cm
#     :return:               a list with parameters (image_siza, lens_focal_distance,
#                            focal_position from transfocator center, divergence of beam after the transfocator)
#     """
#
#     deltas = [(1.0 - xraylib.Refractive_Index_Re(symbol[i],photon_energy_ev*1e-3,density[i])) \
#               for i in range(len(symbol))]
#     focal_f = _transfocator_calculate_focal_distance( deltas=deltas,\
#                                                        nlenses=nlenses_target,radii=nlenses_radii)
#     focal_q = 1.0 / (1.0/focal_f - 1.0/tf_p)
#     div_q = alpha  * lens_diameter / focal_q
#     #corrections
#     source_demagnified = 2.35*sigmaz*focal_q/tf_p
#     if source_demagnified > lens_diameter: source_demagnified = lens_diameter
#     s_target = numpy.sqrt( (div_q*(tf_q-focal_q))**2 + (source_demagnified)**2 )
#     return (s_target,focal_f,focal_q,div_q)
#
#
# def transfocator_nlenses_to_slots(nlenses,nlenses_max=None):
#     """
#     converts the transfocator configuration from a list of the number of lenses of each type,
#     into a list of active (1) or inactive (0) actuators for the slots.
#
#     :param nlenses: the list with number of lenses (e.g., [5,2,0]
#     :param nlenses_max: the maximum number of lenses of each type, usually powers of two minus one.
#                         E.g. [15,3,1]
#     :return: a list of on (1) and off (0) slots, e.g., [1, 0, 1, 0, 0, 1, 0]
#             (first type: 1*1+0*2+1*4+0*8=5, second type: 0*1+1*2=2, third type: 0*1=0)
#     """
#     if nlenses_max == None:
#         nlenses_max = nlenses
#
#     ss = []
#     for i,iopt in enumerate(nlenses):
#         if iopt > nlenses_max[i]:
#             print("Error: i:%d, nlenses: %d, nlenses_max: %d"%(i,iopt,nlenses_max[i]))
#         ncharacters = len("{0:b}".format(nlenses_max[i]))
#
#         si = list( ("{0:0%db}"%(ncharacters)).format(int(iopt)) )
#         si.reverse()
#         ss += si
#
#     on_off = [int(i) for i in ss]
#     #print("transfocator_nlenses_to_slots: nlenses_max: ",nlenses_max," nlenses: ",nlenses," slots: ",on_off)
#     return on_off
#
#
# def _transfocator_calculate_focal_distance(deltas=[0.999998],nlenses=[1],radii=[500e-4]):
#
#     inverse_focal_distance = 0.0
#     for i,nlensesi in enumerate(nlenses):
#         if nlensesi > 0:
#             focal_distance_i = radii[i] / (2.*nlensesi*deltas[i])
#             inverse_focal_distance += 1.0/focal_distance_i
#     if inverse_focal_distance == 0:
#         return 99999999999999999999999999.
#     else:
#         return 1.0/inverse_focal_distance
#
#
# def _tansfocator_guess_focal_position( s_target, p=6500., q=13500.0, sigmaz=6.6e-4, \
#                                       alpha=0.66, lens_diameter=0.05, method=2):
#     x = 1e15
#     if method == 1: # simple sum
#         AA = 2.35*sigmaz/p
#         BB = -(s_target + alpha * lens_diameter)
#         CC = alpha*lens_diameter*q
#
#         cc = numpy.roots([AA,BB,CC])
#         x = cc[1]
#         return x
#
#     if method == 2: # sum in quadrature
#         AA = ( (2.35*sigmaz)**2)/(p**2)
#         BB = 0.0
#         CC = alpha**2 * lens_diameter**2 - s_target**2
#         DD = - 2.0 * alpha**2 * lens_diameter**2 * q
#         EE = alpha**2 * lens_diameter**2 * q**2
#         cc = numpy.roots([AA,BB,CC,DD,EE])
#         for i,cci in enumerate(cc):
#             if numpy.imag(cci) == 0:
#                 return numpy.real(cci)
#
#     return x
#
# def _transfocator_guess_configuration(focal_f_target,deltas=[0.999998],nlenses_max=[15],radii=[500e-4]):
#
#     nn = len(nlenses_max)
#
#     ncombinations = (1+nlenses_max[0]) * (1+nlenses_max[1]) * (1+nlenses_max[2])
#
#     icombinations = 0
#     aa = numpy.zeros((3,ncombinations),dtype=int)
#     bb = numpy.zeros(ncombinations)
#     for i0 in range(1+nlenses_max[0]):
#         for i1 in range(1+nlenses_max[1]):
#             for i2 in range(1+nlenses_max[2]):
#                 aa[0,icombinations] = i0
#                 aa[1,icombinations] = i1
#                 aa[2,icombinations] = i2
#                 bb[icombinations] = focal_f_target - _transfocator_calculate_focal_distance(deltas=deltas,nlenses=[i0,i1,i2],radii=radii)
#                 icombinations += 1
#     bb1 = numpy.abs(bb)
#     ibest = bb1.argmin()
#
#     return (aa[:,ibest]).tolist()
#

#
#
#

def id18_ray_tracing(emittH=4e-9,emittV=1e-11,betaH=35.6,betaV=3.0,number_of_rays=50000,\
                      density=1.845,symbol="Be",tf_p=1000.0,tf_q=1000.0,lens_diameter=0.05,\
                      slots_max=None,slots_on_off=None,photon_energy_ev=14000.0,\
                      slots_lens_thickness=None,slots_steps=None,slots_radii=None,\
                      s_target=10e-4,focal_f=10.0,focal_q=10.0,div_q=1e-6):

    #=======================================================================================================================
    # Gaussian undulator source
    #=======================================================================================================================


    import Shadow
    #import Shadow.ShadowPreprocessorsXraylib as sx

    sigmaXp = numpy.sqrt(emittH/betaH)
    sigmaZp = numpy.sqrt(emittV/betaV)
    sigmaX = emittH/sigmaXp
    sigmaZ = emittV/sigmaZp

    print("\n\nElectron sizes H:%f um, V:%fu m;\nelectron divergences: H:%f urad, V:%f urad"%\
          (sigmaX*1e6, sigmaZ*1e6, sigmaXp*1e6, sigmaZp*1e6))

    # set Gaussian source
    src = Shadow.Source()
    src.set_energy_monochromatic(photon_energy_ev)
    src.set_gauss(sigmaX*1e2,sigmaZ*1e2,sigmaXp,sigmaZp)

    print("\n\nElectron sizes stored H:%f um, V:%f um;\nelectron divergences: H:%f urad, V:%f urad"%\
          (src.SIGMAX*1e4,src.SIGMAZ*1e4,src.SIGDIX*1e6,src.SIGDIZ*1e6))

    src.apply_gaussian_undulator(undulator_length_in_m=2.8, user_unit_to_m=1e-2, verbose=1)

    print("\n\nElectron sizes stored (undulator) H:%f um, V:%f um;\nelectron divergences: H:%f urad, V:%f urad"%\
          (src.SIGMAX*1e4,src.SIGMAZ*1e4,src.SIGDIX*1e6,src.SIGDIZ*1e6))

    print("\n\nSource size in vertical FWHM: %f um\n"%\
          (2.35*src.SIGMAZ*1e4))

    src.NPOINT = number_of_rays
    src.ISTAR1 = 0 # 677543155


    src.write("start.00")

    # create source
    beam = Shadow.Beam()
    beam.genSource(src)
    beam.write("begin.dat")
    src.write("end.00")


    #=======================================================================================================================
    # complete the (detailed) transfocator description
    #=======================================================================================================================


    print("\nSetting detailed Transfocator for ID30B")



    slots_nlenses = numpy.array(slots_max)*numpy.array(slots_on_off)
    slots_empty = (numpy.array(slots_max)-slots_nlenses)

    #
    ####interactive=True, SYMBOL="SiC",DENSITY=3.217,FILE="prerefl.dat",E_MIN=100.0,E_MAX=20000.0,E_STEP=100.0


    Shadow.ShadowPreprocessorsXraylib.prerefl(interactive=False,E_MIN=2000.0,E_MAX=55000.0,E_STEP=100.0,\
                                              DENSITY=density,SYMBOL=symbol,FILE="Be2_55.dat" )

    nslots = len(slots_max)
    prerefl_file = ["Be2_55.dat" for i in range(nslots)]


    print("slots_max:     ",slots_max)
    #print("slots_target:  ",slots_target)
    print("slots_on_off:  ",slots_on_off)
    print("slots_steps:   ",slots_steps)
    print("slots_radii:   ",slots_radii)
    print("slots_nlenses: ",slots_nlenses)
    print("slots_empty:   ",slots_empty)


    #calculate distances, nlenses and slots_empty

    # these are distances p and q with TF length removed
    tf_length = numpy.array(slots_steps).sum()  #tf length in cm
    tf_fs_before = tf_p - 0.5*tf_length     #distance from source to center of transfocator
    tf_fs_after  = tf_q - 0.5*tf_length     # distance from center of transfocator to image

    # for each slot, these are the empty distances before and after the lenses
    tf_p0 = numpy.zeros(nslots)
    tf_q0 = numpy.array(slots_steps) - (numpy.array(slots_max) * slots_lens_thickness)
    # add now the p q distances
    tf_p0[0]  += tf_fs_before
    tf_q0[-1] += tf_fs_after

    print("tf_p0:   ",tf_p0)
    print("tf_q0:   ",tf_q0)
    print("tf_length: %f cm"%(tf_length))


    # build transfocator
    tf = Shadow.CompoundOE(name='TF ID30B')


    tf.append_transfocator(tf_p0.tolist(), tf_q0.tolist(), \
                    nlenses=slots_nlenses.tolist(), radius=slots_radii, slots_empty=slots_empty.tolist(),\
                    thickness=slots_lens_thickness, prerefl_file=prerefl_file,\
                    surface_shape=4, convex_to_the_beam=0, diameter=lens_diameter,\
                    cylinder_angle=0.0,interthickness=50e-4,use_ccc=0)


    itmp = input("SHADOW Source complete. Do you want to run SHADOR trace? [1=Yes,0=No]: ")
    if str(itmp) != "1":
        return

    #trace system
    tf.dump_systemfile()
    beam.traceCompoundOE(tf,write_start_files=0,write_end_files=0,write_star_files=0, write_mirr_files=0)

    #write only last result file
    beam.write("star_tf.dat")
    print("\nFile written to disk: star_tf.dat")



    #
    # #ideal calculations
    #



    print("\n\n\n")
    print("=============================================== TRANSFOCATOR OUTPUTS ==========================================")

    print("\nTHEORETICAL results: ")

    print("REMIND-----With these lenses we obtained (analytically): ")
    print("REMIND-----  focal_f: %f cm"%(focal_f))
    print("REMIND-----  focal_q: %f cm"%(focal_q))
    print("REMIND-----  s_target: %f um"%(s_target*1e4))


    demagnification_factor = tf_p/focal_q
    theoretical_focal_size = src.SIGMAZ*2.35/demagnification_factor



    # analyze shadow results
    print("\nSHADOW results: ")
    st1 = beam.get_standard_deviation(3,ref=0)
    st2 = beam.get_standard_deviation(3,ref=1)

    print("  stDev*2.35: unweighted: %f um, weighted: %f um "%(st1*2.35*1e4,st2*2.35*1e4))

    tk = beam.histo1(3, nbins=75, ref=1, nolost=1, write="HISTO1")

    print("  Histogram FWHM: %f um "%(1e4*tk["fwhm"]))


    print("  Transmitted intensity: %f (source was: %d) (transmission is %f %%) "%(beam.intensity(nolost=1), src.NPOINT, beam.intensity(nolost=1)/src.NPOINT*100))



    #scan around image
    xx1 = numpy.linspace(0.0,1.1*tf_fs_after,11) # position from TF exit plane

    #xx0 = focal_q - tf_length*0.5
    xx0 = focal_q - tf_length*0.5 # position of focus from TF exit plane

    xx2 = numpy.linspace(xx0-100.0,xx0+100,21) # position from TF exit plane
    xx3 = numpy.array([tf_fs_after])

    xx = numpy.concatenate(([-0.5*tf_length],xx1,xx2,[tf_fs_after]))

    xx.sort()


    f = open("id30b.spec","w")
    f.write("#F id30b.spec\n")
    f.write("\n#S 1 calculations for id30b transfocator\n")
    f.write("#N 8\n")
    labels = "  %18s  %18s   %18s   %18s    %18s    %18s    %18s    %18s"%\
             ("pos from source","pos from image","[pos from TF]", "pos from TF center", "pos from focus",\
              "fwhm shadow(stdev)","fwhm shadow(histo)","fwhm theoretical")

    f.write("#L "+labels+"\n")

    out = numpy.zeros((8,xx.size))
    for i,pos in enumerate(xx):
        beam2 = beam.duplicate()
        beam2.retrace(-tf_fs_after+pos)
        fwhm1 = 2.35*1e4*beam2.get_standard_deviation(3,ref=1,nolost=1)
        tk = beam2.histo1(3, nbins=75, ref=1, nolost=1)
        fwhm2 = 1e4*tk["fwhm"]
        #fwhm_th = 1e4*transfocator_calculate_estimated_size(pos,diameter=diameter,focal_distance=focal_q)
        fwhm_th2 = 1e4*numpy.sqrt( (div_q*(pos+0.5*tf_length-focal_q))**2 + theoretical_focal_size**2 )
        #fwhm_th2 = 1e4*( numpy.abs(div_q*(pos-focal_q+0.5*tf_length)) + theoretical_focal_size )


        out[0,i] = tf_fs_before+tf_length+pos
        out[1,i] = -tf_fs_after+pos
        out[2,i] = pos
        out[3,i] = pos+0.5*tf_length
        out[4,i] = pos+0.5*tf_length-focal_q


        out[5,i] = fwhm1
        out[6,i] = fwhm2
        out[7,i] = fwhm_th2


        f.write(" %18.3f   %18.3f   %18.3f   %18.3f    %18.3f    %18.3f    %18.3f    %18.3f \n"%\
                (tf_fs_before+tf_length+pos,\
                 -tf_fs_after+pos,\
                 pos,\
                 pos+0.5*tf_length,\
                 pos+0.5*tf_length-focal_q,\
                 fwhm1,fwhm2,fwhm_th2))

    f.close()
    print("File with beam evolution written to disk: id30b.spec")


    #
    # plots
    #
    itmp = input("Do you want to plot the intensity distribution and beam evolution? [1=yes,0=No]")
    if str(itmp) != "1":
        return

    import matplotlib.pylab as plt

    plt.figure(1)
    plt.plot(out[1,:],out[5,:],'blue',label="fwhm shadow(stdev)")
    plt.plot(out[1,:],out[6,:],'green',label="fwhm shadow(histo1)")
    plt.plot(out[1,:],out[7,:],'red',label="fwhm theoretical")

    plt.xlabel("Distance from image plane [cm]")
    plt.ylabel("spot size [um] ")
    ax = plt.subplot(111)
    ax.legend(bbox_to_anchor=(1.1, 1.05))

    print("Kill graphic to continue.")
    plt.show()

    Shadow.ShadowTools.histo1(beam,3,nbins=75,ref=1,nolost=1,calculate_widths=1)

    input("<Enter> to finish.")
    return None

def id18_full_simulation(photon_energy_ev=14000.0,s_target=20.0e-6,nlenses_target=None):


    if nlenses_target == None:
        force_nlenses = 0
    else:
        force_nlenses = 1

    #
    # define lens setup (general)
    #
    xrl_symbol = ["Be","Be","Be"]
    xrl_density = [1.845,1.845,1.845]
    lens_diameter = 0.05
    nlenses_max = [15,3,1]
    nlenses_radii = [500e-6,1000e-6,1500e-6]

    sigmaz=6.6e-6
    alpha = 0.55

    tf_p = 164.0 # position of the TF measured from the center of the transfocator
    tf_q = 36.0  # position of the image plane measured from the center of the transfocator


    if s_target < 2.35*sigmaz*tf_q/tf_p:
        print("Source size FWHM is: %f um"%(1e6*2.35*sigmaz))
        print("Maximum Demagnifications is: %f um"%(tf_p/tf_q))
        print("Minimum possible size is: %f um"%(1e6*2.35*sigmaz*tf_q/tf_p))
        print("Error: redefine size")
        return


    print("================================== TRANSFOCATOR INPUTS ")
    print("Photon energy: %f eV"%(photon_energy_ev))
    if force_nlenses:
        print("Forced_nlenses: ",nlenses_target)
    else:
        print("target size: %f um"%(1e6*s_target))

    print("materials: ",xrl_symbol)
    print("densities: ",xrl_density)
    print("Lens diameter: %f cm"%(lens_diameter))
    print("nlenses_max:",nlenses_max,"nlenses_radii: ",nlenses_radii)

    print("Source size (sigma): %f um, FWHM: %f um"%(1e6*sigmaz,2.35*1e6*sigmaz))
    print("Distances: tf_p: %f cm, tf_q: %f cm"%(tf_p,tf_q))
    print("alpha: %f"%(alpha))

    print("========================================================")


    if force_nlenses != 1:
        nlenses_target =  transfocator_compute_configuration(photon_energy_ev,s_target,\
            symbol=xrl_symbol,density=xrl_density,\
            nlenses_max=nlenses_max, nlenses_radii=nlenses_radii, lens_diameter=lens_diameter, \
            sigmaz=sigmaz, alpha=alpha, \
            tf_p=tf_p,tf_q=tf_q, verbose=1)


    (s_target,focal_f,focal_q,div_q) = \
        transfocator_compute_parameters(photon_energy_ev, nlenses_target,\
                                        symbol=xrl_symbol,density=xrl_density,\
                                        nlenses_max=nlenses_max, nlenses_radii=nlenses_radii, \
                                        lens_diameter=lens_diameter,\
                                        sigmaz=sigmaz, alpha=alpha,\
                                        tf_p=tf_p,tf_q=tf_q)

    slots_max   = [  1,  2,  4,  8,   1,   2,   1]  # slots
    slots_on_off  = transfocator_nlenses_to_slots(nlenses_target,nlenses_max=nlenses_max)


    print("=============================== TRANSFOCATOR SET")
    #print("deltas: ",deltas)

    if force_nlenses != 1:
        print("nlenses_target (optimized): ",nlenses_target)
    else:
        print("nlenses_target (forced): ",nlenses_target)

    print("With these lenses we obtain: ")
    print("  focal_f: %f m"%(focal_f))
    print("  focal_q: %f m"%(focal_q))
    print("  s_target: %f um"%(s_target*1e6))
    print("  slots_max:    ",slots_max)
    print("  slots_on_off: ",slots_on_off)
    print("==================================================")

    # for theoretical calculations use the focal position and distances given by the target nlenses

    itmp = input("Start SHADOW simulation? [1=yes,0=No]: ")
    if str(itmp) != "1":
        return

    #=======================================================================================================================
    # Inputs
    #=======================================================================================================================


    emittH = 3.9e-9
    emittV = 10e-12
    betaH = 35.6
    betaV = 3.0
    number_of_rays = 50000

    nslots = len(slots_max)
    slots_lens_thickness = [0.3 for i in range(nslots)]   #total thickness of a single lens in cm
    # for each slot, positional gap  of the first lens in cm
    slots_steps    = [  4,   4, 1.9, 6.1,   4,   4, slots_lens_thickness[-1]]
    slots_radii   = [.05, .05, .05, .05, 0.1, 0.1, 0.15]  # radii of the lenses in cm
    AAA= 333
    id18_ray_tracing(emittH=emittH,emittV=emittV,betaH=betaH,betaV=betaV,number_of_rays=number_of_rays,\
                      density=xrl_density[0],symbol=xrl_symbol[0],tf_p=tf_p,tf_q=tf_q,lens_diameter=lens_diameter,\
                      slots_max=slots_max,slots_on_off=slots_on_off,photon_energy_ev=photon_energy_ev,\
                      slots_lens_thickness=slots_lens_thickness,slots_steps=slots_steps,slots_radii=slots_radii,\
                      s_target=s_target,focal_f=focal_f,focal_q=focal_q,div_q=div_q)


def main():

    # this performs the full simulation: calculates the optimum configuration and do the ray-tracing
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

    # itmp = 1 # input("full simulation (ray tracing) ? Enter: \n  0 = No \n  1 = Yes \n> ")
    #
    # if str(itmp) == "1":

    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n")
    id18_full_simulation(photon_energy_ev=photon_energy_kev*1e3,
                          s_target=s_target_um*1e-6,
                          nlenses_target=None)

if __name__ == "__main__":
    main()