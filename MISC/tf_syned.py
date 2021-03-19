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
        SynedObject.__init__(self)
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

    def get_focal_distance_one_lens(self, index, photon_energy=10000.0, verbose=False):
        #
        if verbose: print("\n\n\n ==========  parameters in use : ")

        refraction_index_delta, att_coefficient = \
            self._stack_list[index].get_refraction_index(photon_energy=photon_energy)

        # this is for info...
        number_of_curved_surfaces = self._stack_list[index]._keywords_at_creation["number_of_curved_surfaces"]
        lens_radius = self._stack_list[index]._keywords_at_creation["radius"]
        n_lenses = self._stack_list[index]._keywords_at_creation["n_lenses"]

        if verbose:
            print("\n\nRadius of curvature R = %g um" % (1e6 * lens_radius))
            print("Number of lenses N: %d" % n_lenses)
            print("Number of curved refractive surfaces in a lens Nd = %d" % (number_of_curved_surfaces))
        if number_of_curved_surfaces != 0:
            F = lens_radius / (number_of_curved_surfaces * n_lenses * refraction_index_delta)
            if verbose:
                print("Focal distance F = R / (Nd N delta) = %g m" % (F))
        return F

    def get_focal_distances(self, photon_energy=10000.0): # , mask=None):
        n = self.get_number_of_lenses()
        # if mask is None:
        #     mask = numpy.ones(n)
        out = numpy.zeros(n)

        for i in range(n):
            out[i] = self.get_focal_distance_one_lens(i, photon_energy=photon_energy) #* mask[i]
        return numpy.array(out)

    def get_inverse_focal_distances(self, photon_energy=10000.0): #, mask=None):
        return 1.0 / self.get_focal_distances(photon_energy=photon_energy)

    def get_inverse_focal_distance(self, photon_energy=10000.0, mask=None):
        n = self.get_number_of_lenses()
        if mask is None:
            mask = numpy.ones(n)
        f_inv_out = 0.0
        fs = self.get_inverse_focal_distances(photon_energy=photon_energy)
        for i in range(n):
            if mask[i]: f_inv_out += fs[i]
        return f_inv_out

    def get_focal_distance(self, photon_energy=10000.0, mask=None):
        return 1.0 / self.get_inverse_focal_distance(photon_energy=photon_energy, mask=mask)


    def get_mask_from_index(self, index):
        aa = numpy.binary_repr(index, width=self.get_number_of_lenses())
        return numpy.array(list(aa), dtype=int)

    def get_index_from_mask(self, mask):
        idx = 0
        n = len(mask)
        for i in range(n):
            idx += mask[i] * 2**(n-i-1)
        return idx


    def calculate_focal_length_for_all_masks(self, photon_energy=10000.0):
        f_inverse = self.get_inverse_focal_distances(photon_energy=photon_energy)
        print(">>>>", f_inverse.shape )
        n = self.get_number_of_lenses()
        Fs = numpy.zeros(2**n)
        # print("n=%d, analyzing 2**N=%d configurations: " % (n, 2**n) )
        for i in range(2**n):
            # aa = str(numpy.binary_repr(i, width=n))
            # mask = numpy.array(list(aa), dtype=int)
            mask = self.get_mask_from_index(i)
            Fs[i] = 1.0 / (f_inverse * mask).sum()
            print(i, mask, Fs[i], self.get_index_from_mask(mask))
        return Fs

    def guess_configuration_index(self, p=100.0, q=None, photon_energy=10000.0):
        if q is None:
            f = p
        else:
            f = 1.0 (1.0 / p + 1.0 / q)

        Fn = self.calculate_focal_length_for_all_masks(photon_energy=photon_energy)
        distances = numpy.abs(Fn - f)
        mask_index = numpy.argmin(distances)
        return mask_index


    def guess_configuration_mask(self, p=100.0, q=None, photon_energy=10000.0):
        mask_index = self.guess_configuration_index(p=p, q=q, photon_energy=photon_energy)
        return self.get_mask_from_index(mask_index)

    def info_configuration(self, mask=None):
        n = self.get_number_of_lenses()
        if mask is None:
            mask = numpy.ones(n)
        print("\n\n###############################################################")
        print("AXIS #   OFF(0)/ON(1)   N   MATERIAL   RADIUS")
        for i in range(n):
            lens_radius = self._stack_list[i]._keywords_at_creation["radius"]
            n_lenses = self._stack_list[i]._keywords_at_creation["n_lenses"]
            material = self._stack_list[i].get_material()
            print("%d        %d              %d   %s         %f" % (i, mask[i], n_lenses, material, lens_radius))
        print("###############################################################\n\n")

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
    from srxraylib.plot.gol import plot

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
    l_16x200um = LensBlock(16,  radius=200e-6,   thickness=1e-3)
    l_32x200um = LensBlock(32,  radius=200e-6,   thickness=1e-3)
    l_64x200um = LensBlock(64,  radius=200e-6,   thickness=1e-3)


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

# print(collimating_transfocator.info())
# print(collimating_transfocator.get_focal_distances())
# print(collimating_transfocator.get_focal_distance())
# print(collimating_transfocator.to_dictionary())

Fs = collimating_transfocator.calculate_focal_length_for_all_masks()
print(Fs.min(), Fs.max())

isorted = numpy.argsort(Fs)
x = numpy.arange(Fs.size)
plot(x, Fs[isorted], ylog=True)




idx = collimating_transfocator.guess_configuration_index()
mask = collimating_transfocator.get_mask_from_index(idx)
print(">>>> mask: ", mask.shape, mask)
f_with_mask = collimating_transfocator.get_focal_distance(photon_energy=10000., mask=mask)
print(idx, mask,
        collimating_transfocator.get_focal_distance(photon_energy=10000., mask=mask))

collimating_transfocator.info_configuration(mask=mask)
