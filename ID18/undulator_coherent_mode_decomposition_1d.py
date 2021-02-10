
import numpy as np
from scipy import interpolate


# needed by pySRU
from pySRU.ElectronBeam import ElectronBeam as PysruElectronBeam
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as PysruUndulator
from pySRU.Simulation import create_simulation
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD

def calculate_undulator_emission(
        electron_energy=6.04,
        electron_current=0.2,
        undulator_period=0.032,
        undulator_nperiods=50,
        K=0.25,
        photon_energy=10490.0,
        nsigma=6,
        number_of_points=100,
        distance_to_screen=100,
        scan_direction="V"):


    myelectronbeam = PysruElectronBeam(Electron_energy=electron_energy, I_current=electron_current)
    myundulator = PysruUndulator(K=K, period_length=undulator_period, length=undulator_period * undulator_nperiods)

    theta_central_cone = su.gaussian_central_cone_aperture(ebeam.gamma())
    abscissas = np.linspace(-nsigma * theta_central_cone * distance_to_screen, nsigma * theta_central_cone * distance_to_screen, number_of_points)
    if scan_direction == "H":
        X = abscissas
        Y = np.zeros_like(abscissas)
    elif scan_direction == "V":
        X = np.zeros_like(abscissas)
        Y = abscissas

    print("Calculating energy %g eV" % photon_energy)
    simulation_test = create_simulation(magnetic_structure=myundulator, electron_beam=myelectronbeam,
                                        magnetic_field=None, photon_energy=photon_energy,
                                        traj_method=TRAJECTORY_METHOD_ANALYTIC, Nb_pts_trajectory=None,
                                        rad_method=RADIATION_METHOD_APPROX_FARFIELD, initial_condition=None,
                                        distance=distance_to_screen,
                                        X=X, Y=Y, XY_are_list=True)

    # TODO: this is not nice: I redo the calculations because I need the electric vectors to get polarization
    #       this should be avoided after refactoring pySRU to include electric field in simulations!!
    electric_field = simulation_test.radiation_fact.calculate_electrical_field(
        simulation_test.trajectory, simulation_test.source, X, Y, distance_to_screen)

    E = electric_field._electrical_field
    pol_deg1 = (np.abs(E[:,0]) / (np.abs(E[:,0]) + np.abs(E[:,1]))).flatten() # SHADOW definition!!

    intens1 = simulation_test.radiation.intensity.copy()

    #  Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
    intens1 *= (distance_to_screen * 1e3) ** 2 # photons/mm^2 -> photons/rad^2
    intens1 /= 1e-3 * photon_energy # photons/o.1%bw -> photons/eV

    # unpack trajectory
    T0 = simulation_test.trajectory
    T = np.vstack((T0.t,T0.x,T0.y,T0.z,T0.v_x,T0.v_y,T0.v_z,T0.a_x,T0.a_y,T0.a_z))



    return {'intensity':intens1,
            'polarization':pol_deg1,
            'electric_field':E,
            'trajectory':T,
            'photon_energy': photon_energy,
            "abscissas":abscissas,
            "D":distance_to_screen,
            "theta": abscissas / distance_to_screen,
            }

def backpropagate(input_wavefront,distance=-100.0,magnification_x=1.0):
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

    from srxraylib.plot.gol import plot, plot_image
    plot_from_oe = 100  # set to a large number to avoid plots

    # ##########  SOURCE ##########
    #
    # #
    # # create output_wavefront
    # #
    # #
    # output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
    #                                                                       number_of_points=1000)
    # output_wavefront.set_photon_energy(10000)
    # output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=0, shift=0, beta=0.0922395)
    #
    # if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

    ##########  OPTICAL SYSTEM ##########

    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    # input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 35 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=distance, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x',magnification_x)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    #
    # ---- plots -----
    #
    if plot_from_oe <= 1: plot(output_wavefront.get_abscissas()*1e6, output_wavefront.get_intensity(),
                               title='OPTICAL ELEMENT NR 1',xtitle="x [um]")

    return output_wavefront

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image

    from syned.storage_ring.electron_beam import ElectronBeam
    from syned.storage_ring.magnetic_structures.undulator import Undulator


    # definitions with syned compatibility
    ebeam = ElectronBeam(energy_in_GeV=6.0, current = 0.2)
    su = Undulator.initialize_as_vertical_undulator(K=1.191085, period_length=0.02, periods_number=100)
    photon_energy = su.resonance_energy(ebeam.gamma(),harmonic=1)
    print("Resonance energy: ", photon_energy)
    # other inputs
    distance_to_screen = 100.0
    number_of_points = 1000

    # calculate
    out = calculate_undulator_emission(electron_energy=ebeam.energy(),
                                       electron_current=ebeam.current(),
                                       undulator_period=su.period_length(),
                                       undulator_nperiods=su.number_of_periods(),
                                       K=su.K(),
                                       photon_energy=photon_energy,
                                       nsigma=4,
                                       number_of_points=number_of_points,
                                       distance_to_screen=distance_to_screen,
                                       scan_direction="V")

    # see results
    for key in out.keys():
        print(key)



    # plot(out["abscissas"],out["intensity"],xtitle="x (at %g m)" % distance_to_screen)
    # theta = np.linspace(-30e-6, 30e-6, 200)
    # plot(theta, out["ftheta"](theta),title="interpolated",xtitle="angle")
    #
    # print(out["electric_field"].shape)
    # plot(out["abscissas"], np.abs(out["electric_field"][:, 0]),
    #      out["abscissas"], np.abs(out["electric_field"][:, 1]),
    #      out["abscissas"], np.abs(out["electric_field"][:, 2]),
    #      xtitle="x (at 100m)", ytitle="E", legend=["Ex","Ey","Ez"])

    # backpropagate
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D


    input_wavefront = GenericWavefront1D.initialize_wavefront_from_arrays(out["abscissas"],
                                                                          out["electric_field"][:, 0])
    input_wavefront.set_photon_energy(photon_energy=photon_energy)

    plot(input_wavefront.get_abscissas() * 1e6, input_wavefront.get_intensity(),
         title='SOURCE', xtitle="x [um]")

    output_wavefront = backpropagate(input_wavefront=input_wavefront,
                                     distance=-distance_to_screen,
                                     magnification_x=1.0/distance_to_screen)



    xx = np.linspace(output_wavefront.get_abscissas().min(), output_wavefront.get_abscissas().max(), 200)
    cc1 = output_wavefront.get_interpolated_complex_amplitudes( xx )

    CC1 = np.outer(np.conjugate(cc1), cc1)
    YY = np.outer(np.ones_like(xx), xx)
    print(CC1.shape)

    # yy = output_wavefront.get_interpolated_complex_amplitudes( xx )
    plot_image(np.abs(CC1), title="interpolated")
