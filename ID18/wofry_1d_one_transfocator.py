
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

#
#
#
from syned.beamline.shape import *
from oasys.util.oasys_util import get_fwhm

def get_fwhm(histogram, bins):
    quote = numpy.max(histogram)*0.5
    cursor = numpy.where(histogram >= quote)

    if histogram[cursor].size > 1:
        bin_size    = bins[1]-bins[0]
        fwhm        = bin_size*(cursor[0][-1]-cursor[0][0])
        coordinates = (bins[cursor[0][0]], bins[cursor[0][-1]])
    else:
        fwhm = 0.0
        coordinates = None

    return fwhm, quote, coordinates

#
#
#
class Score():
    def __init__(self, scan_variable_name='x'):
        self.reset()
        self.scan_variable_name = scan_variable_name

    def reset(self):
        self.scan_variable_index = 0
        self.scan_variable_value = []
        self.fwhm = []
        self.intensity_at_center = []
        self.intensity_total = []
        self.intensity_peak = []

    def append(self, wf, scan_variable_value=None):
        fwhm, intensity_total, intensity_at_center, intensity_peak = self.process_wavefront(wf)
        self.fwhm.append(fwhm)
        self.intensity_at_center.append(intensity_at_center)
        self.intensity_total.append(intensity_total)
        self.intensity_peak.append(intensity_peak)
        self.scan_variable_index += 1
        if scan_variable_value is None:
            self.scan_variable_value.append(self.scan_variable_index)
        else:
            self.scan_variable_value.append(scan_variable_value)

    def save(self, filename="tmp.dat"):
        f = open(filename, 'w')
        for i in range(len(self.fwhm)):
            f.write("%g %g %g %g %g\n" % (self.scan_variable_value[i],
                                    1e6*self.fwhm[i],
                                    self.intensity_total[i],
                                    self.intensity_at_center[i],
                                    self.intensity_peak[i]))
        f.close()
        print("File written to disk: %s" % filename)

    def plot(self, title=""):

        x = numpy.array(self.scan_variable_value)


        y = numpy.array(self.intensity_at_center)
        plot(x, y, yrange=[0,1.1*y.max()],
             title=title, ytitle="Intensity at center[a.u.]", xtitle=self.scan_variable_name,
             figsize=(15, 4), show=0)

        # y = numpy.array(self.intensity_total)
        # plot(x, y, yrange=[0,1.1*y.max()],
        #      title=title, ytitle="Beam intensity [a.u.]", xtitle=self.scan_variable_name,
        #      figsize=(15, 4), show=0)

        y = numpy.array(self.fwhm)
        plot(x, y, yrange=[0,1.1*y.max()],
             title=title, ytitle="FWHM [um]", xtitle=self.scan_variable_name,
             figsize=(15, 4), show=1)



    @classmethod
    def process_wavefront(cls, wf):
        I = wf.get_intensity()
        x = wf.get_abscissas()

        fwhm, quote, coordinates = get_fwhm(I, x)
        intensity_at_center = I[I.size // 2]
        intensity_total = I.sum() * (x[1] - x[0])
        intensity_peak = I.max()

        return fwhm, intensity_total, intensity_at_center, intensity_peak




def get_wavefront_intensity_fwhm(wf):
    fwhm, quote, coordinates = get_fwhm(wf.get_intensity(), wf.get_abscissas())
    return fwhm

def get_wavefront_intensity_I0(wf):
    I = wf.get_intensity()
    return I[I.size // 2]

def get_wavefront_intensity_ITOTAL(wf):
    I = wf.get_intensity()
    x = wf.get_abscissas()
    return I.sum() * (x[1] - x[0])


def run_wofry_1d(plot_from=0, mode_x=0):
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012,x_max=0.00012,number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05,amplitude=1,mode_x=mode_x,shift=0,beta=0.0922395)


    if plot_from <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')

    ##########  OPTICAL SYSTEM ##########





    ##########  OPTICAL ELEMENT NUMBER 1 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 35 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=35.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 8.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')


    #
    #---- plots -----
    #
    if plot_from <= 1:  plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 1')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########



    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *

    # boundary_shape=Rectangle(-1.25e-05, 1.25e-05, -1.25e-05, 1.25e-05)
    boundary_shape = Rectangle(-1e-6*slit_size_in_um/2, 1e-6*slit_size_in_um/2, -1e-6*slit_size_in_um/2, 1e-6*slit_size_in_um/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D, WOGaussianSlit1D
    if use_gaussian_slits:
        optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)
    else:
        optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from <= 2: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 2')

    ##########  OPTICAL ELEMENT NUMBER 3 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 30 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=30.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.5)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='INTEGRAL_1D')


    #
    #---- plots -----
    #
    if plot_from <= 3: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 3')

    ##########  OPTICAL ELEMENT NUMBER 4 ##########



    input_wavefront = output_wavefront.duplicate()


    if use_real_lens:
        from orangecontrib.esrf.wofry.util.lens import WOLens1D
        optical_element = WOLens1D.create_from_keywords(
            name='',
            shape=1,
            radius=0.000192435,
            lens_aperture=0.001,
            wall_thickness=5e-05,
            material='Be',
            refraction_index_delta=5.3e-07, # used if material='External'
            att_coefficient=0.00357382, # used if material='External'
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
    else:
        from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D
        optical_element = WOIdealLens1D(name='IdealLensF=28.2',focal_length=28.200000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)


    #
    #---- plots -----
    #
    if plot_from <= 4: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 4')

    ##########  OPTICAL ELEMENT NUMBER 5 ##########



    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_after 99 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=0.000000,    q=q5,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #

    propagation_parameters.set_additional_parameters('magnification_x', magnification_x)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')


    #
    #---- plots -----
    #
    if plot_from <= 5: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 5')

    return output_wavefront

def run_multimode(up_to_mode=0):

    for i in range(up_to_mode+1):
        wf = run_wofry_1d(plot_from=1000, mode_x=i)
        if i == 0:
            WF = wf.duplicate()
        else:
            intens = WF.get_intensity()
            intens += wf.get_intensity()
            WF.set_complex_amplitude = numpy.sqrt(intens)

    return WF

if __name__ == "__main__":


    use_gaussian_slits = False
    use_real_lens = True
    npoints = 100
    do_plot = False
    save_file = True

    NSIGMAS = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 4.0, 6.0]


    for nsigmas in NSIGMAS:
        sc = Score(scan_variable_name='q5 [m]')

        if nsigmas <= 2:
            magnification_x = 2
        if nsigmas <= .5:
            magnification_x = 1
        else:
            magnification_x = 5

        slit_size_in_um = 125.0 / 2.35 * nsigmas
        Q5 = numpy.concatenate((
            numpy.linspace(10, 70, npoints // 3),
            numpy.linspace(71, 300, npoints // 3),  # 49.8098 at source, 470 at slit
            numpy.linspace(301, 500, npoints // 3) ))

        npoints = Q5.size


        if do_plot:
            q5 = Q5[0] + 0.1
            WF1 = run_multimode(up_to_mode=0)
            q5 = 49.8098
            WF2 = run_multimode(up_to_mode=0)
            q5 = 470
            WF3 = run_multimode(up_to_mode=0)
            q5 = Q5[-1]
            WF4 = run_multimode(up_to_mode=0)

            plot(WF1.get_abscissas(), WF1.get_intensity(),
                 WF2.get_abscissas(), WF2.get_intensity(),
                 WF3.get_abscissas(), WF3.get_intensity(),
                 WF4.get_abscissas(), WF4.get_intensity(),
                 legend=["%g"%Q5[0], "%g"%50, "%g"%470, "%g"%Q5[-1]],
                 title="nsigma = %g" % nsigmas)

        for i in range(npoints):
            # run_wofry_1d(plot_from=50)
            q5 = Q5[i]

            WF = run_multimode(up_to_mode=0)
            sc.append(WF, scan_variable_value=q5)

            if numpy.mod(i,10) == 0:
                print(">>>>>> iteration index %d of %d" % (i,npoints))
                # if 0:
                #     plot(WF.get_abscissas(), WF.get_intensity(),
                #          title=">>>>>> iteration index %d of %d" % (i,npoints))

        if save_file:
            if use_real_lens:
                if use_gaussian_slits:
                    sc.save(filename="tmpGR%2.1f.dat" % nsigmas)
                else:
                    sc.save(filename="tmpR%2.1f.dat" % nsigmas)
            else:
                if use_gaussian_slits:
                    sc.save(filename="tmpG%2.1f.dat" % nsigmas)
                else:
                    sc.save(filename="tmp%2.1f.dat" % nsigmas)
        # sc.save(filename="tmpG.dat")
        if do_plot: sc.plot()