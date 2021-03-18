from orangecontrib.esrf.wofry.util.lens import WOLens1D, WOLens
from wofry.beamline.decorators import OpticalElementDecorator
from syned.syned_object import SynedObject
from collections import OrderedDict

import numpy

# lp = dict(aperture=1e-3)
# lp = dict(thickness=1e-3)

def LensBlock(n_lenses=1, radius=500e-6, thickness=5e-5):
    return WOLens1D.create_from_keywords(
        name="Real Lens 1D",
        shape=1,
        radius=radius,
        lens_aperture=0.001,
        wall_thickness=thickness,
        material="Be",  # can be "External"
        refraction_index_delta=5.3e-07,
        att_coefficient=0.00357382,
        number_of_curved_surfaces=2,
        n_lenses=n_lenses,
        error_flag=0,
        error_file="",
        error_edge_management=0,
        write_profile_flag=0,
        write_profile="",
        mis_flag=0,
        xc=0,
        ang_rot=0,
        wt_offset_ffs=0,
        offset_ffs=0,
        tilt_ffs=0,
        wt_offset_bfs=0,
        offset_bfs=0,
        tilt_bfs=0,
    )

class Transfocator(SynedObject, OpticalElementDecorator):
    def __init__(self,
            sr_lensset=None,
            ):
        self._stack_list = []
        if sr_lensset is not None:
            for lens in sr_lensset:
                self._stack_list.append(lens)


    # overwrites the SynedObject method for dealing with list
    def to_dictionary(self):
        dict_to_save = OrderedDict()
        dict_to_save.update({"CLASS_NAME":self.__class__.__name__})

        dict_to_save["multiple_lens_list"] = [el.to_dictionary() for el in self._stack_list]

        return dict_to_save


    def reset(self):
        self._stack_list = []

    def get_number_of_lenses(self):
        return len(self._stack_list)

    def get_focal_distance(self, index, photon_energy=10000.0):
        #
        print("\n\n\n ==========  parameters in use : ")

        refraction_index_delta, att_coefficient = \
            self._stack_list[index].get_refraction_index(photon_energy=photon_energy)

        # this is for info...
        number_of_curved_surfaces = self._stack_list[index]._keywords_at_creation["number_of_curved_surfaces"]
        lens_radius = self._stack_list[index]._keywords_at_creation["radius"]
        n_lenses = self._stack_list[index]._keywords_at_creation["n_lenses"]

        print("\n\nRadius of curvature R = %g um" % (1e6 * lens_radius))
        print("Number of lenses N: %d" % n_lenses)
        print("Number of curved refractive surfaces in a lens Nd = %d" % (number_of_curved_surfaces))
        if number_of_curved_surfaces != 0:
            F = lens_radius / (number_of_curved_surfaces * n_lenses * refraction_index_delta)
            print("Focal distance F = R / (Nd N delta) = %g m" % (F))
        return F

    def get_focal_distances(self, photon_energy=10000.0):
        out = []
        for i in range(self.get_number_of_lenses()):
            out.append(self.get_focal_distance(i))
        return out
    

    # def get_boundaries(self):
    #     boundaries_list = []
    #     for i in range(self.get_number_of_patches()):
    #         boundaries_list.extend(list(self._patch_list[i].get_boundaries()))
    #     return tuple(boundaries_list)
    #
    # def append_patch(self,patch=BoundaryShape()):
    #     self._patch_list.append(patch)

    # def append_rectangle(self,x_left=-0.010,x_right=0.010,y_bottom=-0.020,y_top=0.020):
    #     self.append_patch(Rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top))
    #
    # def append_circle(self,radius, x_center=0.0, y_center=0.0):
    #     self.append_patch(Circle(radius, x_center=x_center, y_center=y_center))
    #
    # def append_ellipse(self,a_axis_min, a_axis_max, b_axis_min, b_axis_max):
    #     self.append_patch(Ellipse(a_axis_min, a_axis_max, b_axis_min, b_axis_max))
    #
    # def get_patches(self):
    #     return self._patch_list
    #
    # def get_patch(self,index):
    #     return self.get_patches()[index]
    #
    # def get_name_of_patch(self,index):
    #     return self._patch_list[index].__class__.__name__


if __name__ == "__main__":
    l_1x10mm =   LensBlock(1,   radius=10000e-6, thickness=1e-3)
    l_1x5mm =    LensBlock(1,   radius=5000e-6,  thickness=1e-3)
    l_1x2mm =    LensBlock(1,   radius=2000e-6,  thickness=1e-3)
    l_1x1mm =    LensBlock(1,   radius=1000e-6,  thickness=1e-3)
    l_1x500um =  LensBlock(1,   radius=500e-6,   thickness=1e-3)
    l_1x500um =  LensBlock(1,   radius=500e-6,   thickness=1e-3)
    l_1x300um =  LensBlock(1,   radius=300e-6,   thickness=1e-3)
    l_1x200um =  LensBlock(1,   radius=200e-6,   thickness=1e-3)
    l_2x200um =  LensBlock(2,   radius=200e-6,   thickness=1e-3)
    l_4x200um =  LensBlock(4,   radius=200e-6,   thickness=1e-3)
    l_8x200um =  LensBlock(8,   radius=200e-6,   thickness=1e-3)
    l_16x200um = LensBlock(16, radius=200e-6,   thickness=1e-3)
    l_32x200um = LensBlock(32, radius=200e-6,   thickness=1e-3)
    l_64x200um = LensBlock(64, radius=200e-6,   thickness=1e-3)


collimating_transfocator = Transfocator(
    (
        l_1x5mm,
        l_1x2mm,
        l_1x1mm,
        l_1x500um,
        l_1x300um,
        l_1x200um,
        l_2x200um,
        l_4x200um,
        l_8x200um,
    )
)

print(collimating_transfocator.info())
print(collimating_transfocator.get_focal_distances())


