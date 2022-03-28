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


#
# SOURCE========================
#


def run_source_h(my_mode_index=0):
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(7000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.05281e-05, amplitude=1, mode_x=0, shift=0, beta=0.127769)
    # previous command is useless but...
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.05281e-05, amplitude=1, mode_x=my_mode_index, shift=0,
                                               beta=0.127769)
    return output_wavefront


#
# SOURCE========================
#


def run_source_v(my_mode_index=0):
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-5e-05, x_max=5e-05,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(7000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=5.84299e-06, amplitude=1, mode_x=0, shift=0, beta=1.56094)
    # previous command is useless but...
    output_wavefront.set_gaussian_hermite_mode(sigma_x=5.84299e-06, amplitude=1, mode_x=my_mode_index, shift=0,
                                               beta=1.56094)
    return output_wavefront

#
# BEAMLINE========================
#


def run_beamline_h(output_wavefront,slit=50e-6, gaussian_slit=True):
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
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=35.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 4.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-slit/2, slit/2, -slit/2, slit/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOGaussianSlit1D, WOSlit1D
    if gaussian_slit:
        optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)
    else:
        optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    return output_wavefront






#
# BEAMLINE========================
#


def run_beamline_v(output_wavefront,slit=50e-6, gaussian_slit=True):
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
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=35.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 10.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-slit/2, slit/2, -slit/2, slit/2)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOGaussianSlit1D, WOSlit1D
    if gaussian_slit:
        optical_element = WOGaussianSlit1D(boundary_shape=boundary_shape)
    else:
        optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    return output_wavefront


#
# MAIN FUNCTION========================
#


def main_h(slit=50e-6,gaussian_slit=False):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(50):
        output_wavefront = run_source_h(my_mode_index=my_mode_index)
        output_wavefront = run_beamline_h(output_wavefront, slit=slit, gaussian_slit=gaussian_slit)
        tally.append(output_wavefront)

    return tally
    # tally.plot_cross_spectral_density()
    # tally.plot_spectral_density()
    # tally.plot_occupation()



def main_v(slit=50e-6,gaussian_slit=False):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(50):
        output_wavefront = run_source_v(my_mode_index=my_mode_index)
        output_wavefront = run_beamline_v(output_wavefront, slit=slit, gaussian_slit=gaussian_slit)
        tally.append(output_wavefront)

    return tally
    # tally.plot_cross_spectral_density()
    # tally.plot_spectral_density()
    # tally.plot_occupation()

# see https://stackoverflow.com/questions/33159134/matplotlib-y-axis-label-with-multiple-colors
def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw))
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.35, -0.13) , #(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw))
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.12, 0.2), #(-0.10, 0.2),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)
#
# MAIN========================
#
if __name__ == "__main__":

    do_calculate = False
    gaussian_slit = False



    #
    #
    #




    if do_calculate:
        if gaussian_slit:
            outfile = "coherent_fraction_vs_slit_apertureG.dat"
        else:
            outfile = "coherent_fraction_vs_slit_aperture.dat"


        # slits = numpy.concatenate((numpy.linspace(10e-6,310e-6, 101), numpy.linspace(320e-6, 0.0015, 21)))
        # slits = numpy.linspace(0.0015, 0.0025, 5)
        slits = numpy.concatenate((numpy.linspace(10e-6,310e-6, 101), numpy.linspace(320e-6, 0.0015, 21), numpy.linspace(0.0015, 0.0025, 5)))


        f = open(outfile, "w")

        for i,slit in enumerate(slits):
            print(">>>>>>>>>>>%d of %d slit: %g" % (i,slits.size-1, slit*1e6))
            tally_v = main_v(slit=slit, gaussian_slit=gaussian_slit)
            tally_h = main_h(slit=slit, gaussian_slit=gaussian_slit)


            modes_v, occ_v = tally_v.get_occupation()
            modes_h, occ_h = tally_h.get_occupation()

            print(">>>>", tally_v.get_cross_pectral_density().shape, tally_v.get_abscissas().shape)
            sd_h = numpy.trapz( tally_v.get_spectral_density(), tally_v.get_abscissas())
            sd_v = numpy.trapz( tally_h.get_spectral_density(), tally_h.get_abscissas())

            print("slit, CF H, V: ", 1e6*slit, occ_h[0], occ_v[0], sd_h, sd_v)
            print("slit, INT H, V: ", 1e6*slit, sd_h, sd_v)
            f.write("%g %g %g %g %g\n" % (slit, occ_h[0], occ_v[0], sd_h, sd_v))
        f.close()
        print("File written to disk: %s" % outfile)


    #
    # plots
    #
    a = numpy.loadtxt("coherent_fraction_vs_slit_apertureG.dat")
    b = numpy.loadtxt("coherent_fraction_vs_slit_aperture.dat")

    from srxraylib.plot.gol import plot, plot_show

    if False:
        g = plot(
             1e6 * a[:,0], a[:,1],
             1e6 * a[:,0], a[:,2],
             1e6 * b[:, 0], b[:, 1],
             1e6 * b[:, 0], b[:, 2],
             legend=['Horizontal Gaussian', 'Vertical Gaussian', 'Horizontal Rectangular', 'Vertical Rectangular'],
             xtitle="Slit aperture [um]", ytitle="Coherent Fraction",
             color = ['green','blue', 'green', 'blue'],
             linestyle=['--','--',None,None],
             xlog=True, yrange=[0,1.01], show=False)

        g[1].grid()
        import matplotlib.pylab as plt
        # locs, labels = plt.yticks()
        plt.yticks(numpy.arange(0, 1.1, step=0.1))
        plot_show()


    #
    # only square aperture
    #

    if True:
        import matplotlib.pylab as plt
        import matplotlib
        matplotlib.rc('xtick', labelsize=14)
        matplotlib.rc('ytick', labelsize=12)
        g = plot(
             1e6 * b[:, 0], b[:, 1],
             1e6 * b[:, 0], b[:, 2],
            1e6 * b[:, 0], b[:, 3] / b[:, 3].max(),
            1e6 * b[:, 0], b[:, 4] / b[:, 4].max(),
             legend=['CF$_x$', 'CF$_y$','Integrated $\mathcal{I}_x$', 'Integrated $\mathcal{I}_y$'],
             xtitle="", # "Slit aperture [$\mu$m]",
             # ytitle="Coherent Fraction or Normalized Integrated Intensity",
             ytitle="",
             color = ['blue','blue','green','green',],
             linestyle=[None,":",None,":"],
             xlog=True, xrange=[0.9e1, 1.1e3], yrange=[0,1.01], show=False)

        # plt.ylabel("Coherent Fraction or Normalized Integrated Intensity")
        # multicolor_ylabel(g[1], ('Line1', 'and', 'Line2', 'with', 'extra', 'colors!'),
        #                   ('r', 'k', 'b', 'k', 'm', 'g'),
        #                   axis='y', size=15, weight='bold')
        multicolor_ylabel(g[1], ('CF', 'and', 'Integrated Intensity'),
                          ('g', 'k', 'b'),
                          axis='y', size=None, weight='bold')
        multicolor_ylabel(g[1], ("Slit aperture [$\mu$m]", ""),
                          ('k','b'),
                          axis='x', size=None, weight='bold')
        g[1].grid()

        # locs, labels = plt.yticks()
        plt.savefig("CFvsGap.pdf")
        plt.yticks(numpy.arange(0, 1.1, step=0.1))
        plot_show()
        print("File written to disk: CFvsGap.pdf")