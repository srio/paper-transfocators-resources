
#
# Import section
#
import numpy

from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

from wofryimpl.propagator.propagators1D.fresnel import Fresnel1D
from wofryimpl.propagator.propagators1D.fresnel_convolution import FresnelConvolution1D
from wofryimpl.propagator.propagators1D.fraunhofer import Fraunhofer1D
from wofryimpl.propagator.propagators1D.integral import Integral1D
from wofryimpl.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
from wofryimpl.propagator.propagators1D.fresnel_zoom_scaling_theorem import FresnelZoomScaling1D

import os

#
# SOURCE========================
#




def run_source(my_mode_index=0):
    global coherent_mode_decomposition
    try:
        if my_mode_index == 0: raise Exception()
        tmp = coherent_mode_decomposition
    except:        
        
        
        ##########  SOURCE ##########
        
        
        #
        # create output_wavefront
        #
        #
        from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import UndulatorCoherentModeDecomposition1D
        coherent_mode_decomposition = UndulatorCoherentModeDecomposition1D(
            electron_energy=6,
            electron_current=0.2,
            undulator_period=0.018,
            undulator_nperiods=138,
            K=1.85108,
            photon_energy=7000,
            abscissas_interval=0.00025,
            number_of_points=1000,
            distance_to_screen=100,
            scan_direction='V',
            sigmaxx=2.97321e-05,
            sigmaxpxp=4.37237e-06,
            useGSMapproximation=False,)
        # make calculation
        coherent_mode_decomposition_results = coherent_mode_decomposition.calculate()
        
        mode_index = 0
        output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(mode_index)
    output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(my_mode_index)
    return output_wavefront


#
# BEAMLINE========================
#




def run_beamline(output_wavefront, radius1=0.001):
    
    
    ##########  OPTICAL SYSTEM ##########
    
    
    
    
    
    ##########  OPTICAL ELEMENT NUMBER 1 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
    
    optical_element = WOScreen1D()
    
    # drift_before 36 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=36.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 5.0)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='INTEGRAL_1D')
    
    
    ##########  OPTICAL ELEMENT NUMBER 2 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape=Rectangle(-2e-05, 2e-05, -2e-05, 2e-05)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)
    
    # no drift in this element 
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    
    
    ##########  OPTICAL ELEMENT NUMBER 3 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
    
    optical_element = WOScreen1D()
    
    # drift_before 29 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=29.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 2.5)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')
    
    
    ##########  OPTICAL ELEMENT NUMBER 4 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.lens import WOLens1D
    
    optical_element = WOLens1D.create_from_keywords(
        name='',
        shape=1,
        radius=radius1, #0.0004125,
        lens_aperture=0.001,
        wall_thickness=5e-05,
        material='Be',
        number_of_curved_surfaces=2,
        n_lenses=1,
        error_flag=0,
        error_file='<none>',
        error_edge_management=0,
        write_profile_flag=0,
        write_profile='profile1D.dat',
        mis_flag=0,
        xc=0,
        ang_rot=0,
        wt_offset_ffs=0,
        offset_ffs=0,
        tilt_ffs=0,
        wt_offset_bfs=0,
        offset_bfs=0,
        tilt_bfs=0)
    
    # no drift in this element 
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    
    
    ##########  OPTICAL ELEMENT NUMBER 5 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
    
    optical_element = WOScreen1D()
    
    # drift_before 127 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=127.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.25)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')




    ###################################  CORRECTOR TO COMPUTE BEST FOCUS #####################################

    if my_mode_index == 0:
        input_wavefront = output_wavefront.duplicate()
        from orangecontrib.esrf.wofry.util.thin_object_corrector import WOThinObjectCorrector1D  # TODO update

        optical_element = WOThinObjectCorrector1D(
            name='',
            file_with_thickness_mesh_flag=0,
            file_with_thickness_mesh='profile1D.dat',
            material='Be',
            focus_at=8,
            wall_thickness=5e-05,
            apply_correction_to_wavefront=0,
            fit_fraction_in_length=0.1,
            fit_filename='tmp.txt')

        print("\n\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", i, my_mode_index)
        os.system("cat tmp.txt")
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n\n")

        # no drift in this element
        output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    ##########################################################################################################
    
    ##########  OPTICAL ELEMENT NUMBER 6 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from orangecontrib.esrf.wofry.util.lens import WOLens1D
    a = numpy.loadtxt("tmp.txt",skiprows=3)
    radius2 = a[0]
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> using radius1, radius2", my_mode_index, radius1, radius2)
    optical_element = WOLens1D.create_from_keywords(
        name='',
        shape=1,
        radius=radius2, #0.0001048,
        lens_aperture=0.001,
        wall_thickness=5e-05,
        material='Be',
        number_of_curved_surfaces=2,
        n_lenses=1,
        error_flag=0,
        error_file='<none>',
        error_edge_management=0,
        write_profile_flag=0,
        write_profile='profile1D.dat',
        mis_flag=0,
        xc=0,
        ang_rot=0,
        wt_offset_ffs=0,
        offset_ffs=0,
        tilt_ffs=0,
        wt_offset_bfs=0,
        offset_bfs=0,
        tilt_bfs=0)
    
    # no drift in this element 
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    
    
    ##########  OPTICAL ELEMENT NUMBER 7 ##########
    
    
    
    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
    
    optical_element = WOScreen1D()
    
    # drift_before 8 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=8.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.1)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='INTEGRAL_1D')
    return output_wavefront


#
# MAIN FUNCTION========================
#



if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes
    import xraylib

    fileroot = "case7keV_8h"
    subdirectory = "./" #trajectories_8h"

    F1 = numpy.linspace(5,100,200) #numpy.array([29.64]) # numpy.linspace(15,60,3) #
    R1 = []
    for F in F1:
        xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", 7,1.85)).real
        R = F * (2 * xrl_delta)
        R1.append(R)
        print("F: %g  R_Be [m]= %g" % (F, R))

    for i in range(len(R1)):
        tally = TallyCoherentModes()
        for my_mode_index in range(50):
            output_wavefront = run_source(my_mode_index=my_mode_index)
            output_wavefront = run_beamline(output_wavefront, radius1=R1[i])
            tally.append(output_wavefront)

        # tally.plot_cross_spectral_density(show=0,filename="%s_cross_spectral_density.png" % fileroot)
        tally.plot_spectral_density(show=0,filename="%s/%s_spectral_density_%03d.png" % (subdirectory, fileroot, i), title="R1: %g" % F1[i])
        # tally.plot_occupation(show=0,filename="%s_occupation.png" % fileroot)

        tally.save_spectral_density(filename="%s/%s_spectral_density_%03d.dat" % (subdirectory, fileroot,i))
        os.system("echo %03i >> tmp.txt" % i)
        os.system("echo %g >> tmp.txt" % R1[i])
        os.system("echo %g >> tmp.txt" % F1[i])
        os.system("mv tmp.txt %s/%s_%03d.txt" % (subdirectory, fileroot, i) )

        tally.save_occupation(filename="%s/%s_occupation_%03d.dat" % (subdirectory, fileroot, i))
