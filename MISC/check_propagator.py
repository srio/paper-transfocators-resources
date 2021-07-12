
#
# Import section
#
import numpy

from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D

from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
from wofryimpl.propagator.propagators2D.fresnel import Fresnel2D
from wofryimpl.propagator.propagators2D.fresnel_convolution import FresnelConvolution2D
from wofryimpl.propagator.propagators2D.fraunhofer import Fraunhofer2D
from wofryimpl.propagator.propagators2D.integral import Integral2D
from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D

from srxraylib.plot.gol import plot, plot_image
plot_from_oe = 10 # set to a large number to avoid plots
import matplotlib.pylab as plt



##########  SOURCE ##########


#
# create output_wavefront
#
#
output_wavefront = GenericWavefront2D.initialize_wavefront_from_range(x_min=-0.00012,x_max=0.00012,
                                                                      y_min=-5e-05,y_max=5e-05,
                                                                      number_of_points=(401,201))
output_wavefront.set_photon_energy(7000)
output_wavefront.set_gaussian_hermite_mode(sigma_x=3.00818e-05,sigma_y=6.99408e-06,amplitude=1,nx=4,ny=1,betax=0.129748,betay=1.01172)


if plot_from_oe <= 0: plot_image(output_wavefront.get_intensity(),output_wavefront.get_coordinate_x(),output_wavefront.get_coordinate_y(),aspect='auto',title='SOURCE')


##########  OPTICAL SYSTEM ##########





##########  OPTICAL ELEMENT NUMBER 1 ##########



input_wavefront = output_wavefront.duplicate()
from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen

optical_element = WOScreen()

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
propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
propagation_parameters.set_additional_parameters('magnification_x', 8.0)
propagation_parameters.set_additional_parameters('magnification_y', 10.0)
#
propagator = PropagationManager.Instance()
try:
    propagator.add_propagator(FresnelZoomXY2D())
except:
    pass
output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_XY_2D')


#
#---- plots -----
#
if plot_from_oe <= 1: plot_image(output_wavefront.get_intensity(),output_wavefront.get_coordinate_x(),output_wavefront.get_coordinate_y(),aspect='auto',title='OPTICAL ELEMENT NR 1')


##########  OPTICAL ELEMENT NUMBER 2 ##########



input_wavefront = output_wavefront.duplicate()
from syned.beamline.shape import Rectangle
boundary_shape=Rectangle(-0.5, 0.5, -0.5, 0.5)
from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit
optical_element = WOSlit(boundary_shape=boundary_shape)

# no drift in this element
output_wavefront = optical_element.applyOpticalElement(input_wavefront)


#
#---- plots -----
#
if plot_from_oe <= 2: plot_image(output_wavefront.get_intensity(),output_wavefront.get_coordinate_x(),output_wavefront.get_coordinate_y(),aspect='auto',title='OPTICAL ELEMENT NR 2')


##########  OPTICAL ELEMENT NUMBER 3 ##########



input_wavefront = output_wavefront.duplicate()
from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen

optical_element = WOScreen()

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
propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
propagation_parameters.set_additional_parameters('magnification_x', 1.5)
propagation_parameters.set_additional_parameters('magnification_y', 2.0)
#
propagator = PropagationManager.Instance()
try:
    propagator.add_propagator(FresnelZoomXY2D())
except:
    pass
output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_XY_2D')


#
#---- plots -----
#
if plot_from_oe <= 3: plot_image(output_wavefront.get_intensity(),output_wavefront.get_coordinate_x(),output_wavefront.get_coordinate_y(),aspect='auto',title='OPTICAL ELEMENT NR 3')



intensity = output_wavefront.get_intensity()
x = output_wavefront.get_coordinate_x()
y = output_wavefront.get_coordinate_y()
i_x = intensity.sum(axis=1)
i_x /= i_x.max()
i_y = intensity.sum(axis=0)
i_y /= i_y.max()

f1, ax2 = plot(x, i_x,
     y, i_y,
     show=0)
ax2.xaxis.grid()
ax2.yaxis.grid()
plt.show()