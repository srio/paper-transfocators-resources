from orangecontrib.esrf.wofry.util.lens import WOLens1D, WOLens
from wofry.beamline.decorators import OpticalElementDecorator
from syned.syned_object import SynedObject
from collections import OrderedDict

import numpy

# lp = dict(aperture=1e-3)
# lp = dict(thickness=1e-3)



class WOTransfocator(SynedObject, OpticalElementDecorator):
    # note that a lens block is a CRL. And it is represented by a WOLens object with n_lenses
    # Therefore the CRL (represented by a lens) has zero interlens-distance (no propagation)
    def __init__(self,
                 lens_block_list=None,
                 lens_block_status=None,
                 ):
        SynedObject.__init__(self)
        self._lens_block_list = []
        if lens_block_list is not None:
            for lens_block in lens_block_list:
                self._lens_block_list.append(lens_block)

        self._lens_block_status = lens_block_status

    # overwrites the SynedObject method for dealing with list
    def to_dictionary(self):
        dict_to_save = OrderedDict()
        dict_to_save.update({"CLASS_NAME":self.__class__.__name__})
        dict_to_save["transfocator"] = [el.to_dictionary() for el in self._lens_block_list]
        return dict_to_save

    def reset(self):
        self._lens_block_list = []
        self._lens_block_status = None

    def get_number_of_lens_blocks(self):
        return len(self._lens_block_list)

    def set_lens_block_status(self, lens_block_status):
        if lens_block_status is None:
            self._lens_block_status = lens_block_status
            return

        try:
            n = lens_block_status.size
        except:
            n = len(lens_block_status)

        if n != self.get_number_of_lens_blocks():
            raise Exception("Bad lens_block_status dimension")
        else:
            self._lens_block_status = numpy.array(lens_block_status)

    def get_lens_block_status(self):
        if self._lens_block_status is None:
            return numpy.ones(self.get_number_of_lens_blocks())
        else:
            return self._lens_block_status

    def append_lens_block(self,
                          n_lenses=1,
                          radius=500e-6,
                          wall_thickness=5e-5,
                          thickness=1e-3,
                          name="Real Lens 1D"):
        self._lens_block_list.append( WOLens1D.create_from_keywords(
            name=name,
            shape=1,
            radius=radius,
            lens_aperture=0.001,
            wall_thickness=wall_thickness,
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
        ))


    def get_focal_distance_one_lens_block(self, index, photon_energy=10000.0, verbose=False):
        #
        if verbose: print("\n\n\n ==========  parameters in use : ")

        refraction_index_delta, att_coefficient = \
            self._lens_block_list[index].get_refraction_index(photon_energy=photon_energy)

        # this is for info...
        number_of_curved_surfaces = self._lens_block_list[index]._keywords_at_creation["number_of_curved_surfaces"]
        lens_radius = self._lens_block_list[index]._keywords_at_creation["radius"]
        n_lenses = self._lens_block_list[index]._keywords_at_creation["n_lenses"]

        if verbose:
            print("\n\nRadius of curvature R = %g um" % (1e6 * lens_radius))
            print("Number of lenses N: %d" % n_lenses)
            print("Number of curved refractive surfaces in a lens Nd = %d" % (number_of_curved_surfaces))
        if number_of_curved_surfaces != 0:
            F = lens_radius / (number_of_curved_surfaces * n_lenses * refraction_index_delta)
            if verbose:
                print("Focal distance F = R / (Nd N delta) = %g m" % (F))
        return F

    def get_focal_distances_of_lens_blocks(self, photon_energy=10000.0):
        n = self.get_number_of_lens_blocks()
        out = numpy.zeros(n)

        for i in range(n):
            out[i] = self.get_focal_distance_one_lens_block(i, photon_energy=photon_energy)
        return numpy.array(out)

    def get_inverse_focal_distances_of_lens_blocks(self, photon_energy=10000.0):
        return 1.0 / self.get_focal_distances_of_lens_blocks(photon_energy=photon_energy)

    def get_inverse_focal_distance(self, photon_energy=10000.0):
        n = self.get_number_of_lens_blocks()
        lens_block_status = self.get_lens_block_status()
        f_inv_out = 0.0
        fs = self.get_inverse_focal_distances_of_lens_blocks(photon_energy=photon_energy)
        for i in range(n):
            if lens_block_status[i]: f_inv_out += fs[i]
        return f_inv_out

    def get_focal_distance(self, photon_energy=10000.0):
        inv_f = self.get_inverse_focal_distance(photon_energy=photon_energy)
        if inv_f == 0:
            return numpy.inf
        else:
            return 1.0 / inv_f


    def get_mask_from_index(self, index):
        aa = numpy.binary_repr(index, width=self.get_number_of_lens_blocks())
        return numpy.array(list(aa), dtype=int)

    def get_index_from_mask(self, mask):
        idx = 0
        n = len(mask)
        for i in range(n):
            idx += mask[i] * 2**(n-i-1)
        return idx

    def get_index_from_current_lens_block_status(self):
        return self.get_index_from_mask(self.get_lens_block_status())


    def calculate_focal_length_for_all_lens_block_status(self, photon_energy=10000.0):
        f_inverse = self.get_inverse_focal_distances_of_lens_blocks(photon_energy=photon_energy)
        n = self.get_number_of_lens_blocks()
        Fs = numpy.zeros(2**n)
        for i in range(2**n):
            mask = self.get_mask_from_index(i)
            Fs[i] = 1.0 / (f_inverse * mask).sum()
            # print(i, mask, Fs[i], self.get_index_from_mask(mask))
        return Fs

    def guess_configuration_index(self, p=100.0, q=None, photon_energy=10000.0):
        if q is None:
            f = p
        else:
            f = 1.0 (1.0 / p + 1.0 / q)
        Fn = self.calculate_focal_length_for_all_lens_block_status(photon_energy=photon_energy)
        distances = numpy.abs(Fn - f)
        mask_index = numpy.argmin(distances)
        return mask_index

    def guess_lens_block_status(self, p=100.0, q=None, photon_energy=10000.0):
        mask_index = self.guess_configuration_index(p=p, q=q, photon_energy=photon_energy)
        return self.get_mask_from_index(mask_index)

    def guess_and_set_lens_block_status(self, p=100.0, q=None, photon_energy=10000.0):
        self.set_lens_block_status(self.guess_lens_block_status(p=p, q=q, photon_energy=photon_energy))


    def info(self):
        txt = ""
        for lens in self._lens_block_list:
            txt += lens.info()
        return txt

    def info_configuration(self, photon_energy=None):
        n = self.get_number_of_lens_blocks()
        lens_block_status = self.get_lens_block_status()

        txt = "\n\n###############################################################"
        txt += "\nBLOCK    STATUS                                  "
        txt += "\nAXIS #   OFF(0)/ON(1)   N   MATERIAL   RADIUS [m]"
        for i in range(n):
            lens_radius = self._lens_block_list[i]._keywords_at_creation["radius"]
            n_lenses = self._lens_block_list[i]._keywords_at_creation["n_lenses"]
            material = self._lens_block_list[i].get_material()
            txt += "\n%d        %d              %d   %s         %f" % (i, lens_block_status[i], n_lenses, material, lens_radius)
        if photon_energy is not None:
            txt += "\nFocal distance [m] at E=%g eV) is: " % \
                  (photon_energy) + repr(self.get_focal_distance(photon_energy))
        txt += "\n###############################################################\n\n"

        return txt


if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    # l_1x10mm =   LensBlock(1,   radius=10000e-6, thickness=1e-3, name="Transfocator crl l_1x10mm")
    # l_1x5mm =    LensBlock(1,   radius=5000e-6,  thickness=1e-3, name="Transfocator crl l_1x5mm")
    # l_1x2mm =    LensBlock(1,   radius=2000e-6,  thickness=1e-3, name="Transfocator crl l_1x2mm")
    # l_1x1mm =    LensBlock(1,   radius=1000e-6,  thickness=1e-3, name="Transfocator crl l_1x1mm")
    # l_1x500um =  LensBlock(1,   radius=500e-6,   thickness=1e-3, name="Transfocator crl l_1x500um")
    # l_1x500um =  LensBlock(1,   radius=500e-6,   thickness=1e-3, name="Transfocator crl l_1x500um")
    # l_1x300um =  LensBlock(1,   radius=300e-6,   thickness=1e-3, name="Transfocator crl l_1x300um")
    # l_1x200um =  LensBlock(1,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_1x200um")
    # l_2x200um =  LensBlock(2,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_2x200um")
    # l_4x200um =  LensBlock(4,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_4x200um")
    # l_8x200um =  LensBlock(8,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_8x200um")
    # l_16x200um = LensBlock(16,  radius=200e-6,   thickness=1e-3, name="Transfocator crl l_16x200um")
    # l_32x200um = LensBlock(32,  radius=200e-6,   thickness=1e-3, name="Transfocator crl l_32x200um")
    # l_64x200um = LensBlock(64,  radius=200e-6,   thickness=1e-3, name="Transfocator crl l_64x200um")

    l_1x10mm  =  dict(n_lenses=1,   radius=10000e-6, thickness=1e-3, name="Transfocator crl l_1x10mm")
    l_1x5mm   =  dict(n_lenses=1,   radius=5000e-6,  thickness=1e-3, name="Transfocator crl l_1x5mm")
    l_1x2mm   =  dict(n_lenses=1,   radius=2000e-6,  thickness=1e-3, name="Transfocator crl l_1x2mm")
    l_1x1mm   =  dict(n_lenses=1,   radius=1000e-6,  thickness=1e-3, name="Transfocator crl l_1x1mm")
    l_1x500um =  dict(n_lenses=1,   radius=500e-6,   thickness=1e-3, name="Transfocator crl l_1x500um")
    l_1x500um =  dict(n_lenses=1,   radius=500e-6,   thickness=1e-3, name="Transfocator crl l_1x500um")
    l_1x300um =  dict(n_lenses=1,   radius=300e-6,   thickness=1e-3, name="Transfocator crl l_1x300um")
    l_1x200um =  dict(n_lenses=1,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_1x200um")
    l_2x200um =  dict(n_lenses=2,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_2x200um")
    l_4x200um =  dict(n_lenses=4,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_4x200um")
    l_8x200um =  dict(n_lenses=8,   radius=200e-6,   thickness=1e-3, name="Transfocator crl l_8x200um")
    l_16x200um = dict(n_lenses=16,  radius=200e-6,   thickness=1e-3, name="Transfocator crl l_16x200um")
    l_32x200um = dict(n_lenses=32,  radius=200e-6,   thickness=1e-3, name="Transfocator crl l_32x200um")
    l_64x200um = dict(n_lenses=64,  radius=200e-6,   thickness=1e-3, name="Transfocator crl l_64x200um")

collimating_transfocator = WOTransfocator()

# collimating_transfocator.append_lens_block(1,radius=5000e-6,  thickness=1e-3, name="Transfocator crl l_1x5mm")
# collimating_transfocator.append_lens_block(1,radius=2000e-6,  thickness=1e-3, name="Transfocator crl l_1x2mm")
# collimating_transfocator.append_lens_block(1,radius=1000e-6,  thickness=1e-3, name="Transfocator crl l_1x1mm")
# collimating_transfocator.append_lens_block(1,radius=500e-6,  thickness=1e-3, name="Transfocator crl  l_1x500um")
# collimating_transfocator.append_lens_block(1,radius=300e-6,  thickness=1e-3, name="Transfocator crl  l_1x300um")
# collimating_transfocator.append_lens_block(1,radius=200e-6,  thickness=1e-3, name="Transfocator crl  l_1x200um")
# collimating_transfocator.append_lens_block(2,radius=200e-6,  thickness=1e-3, name="Transfocator crl  l_2x200um")
# collimating_transfocator.append_lens_block(4,radius=200e-6,  thickness=1e-3, name="Transfocator crl  l_4x200um")
# collimating_transfocator.append_lens_block(8,radius=200e-6,  thickness=1e-3, name="Transfocator crl  l_8x200um")


collimating_transfocator.append_lens_block(**l_1x10mm  )
collimating_transfocator.append_lens_block(**l_1x5mm   )
collimating_transfocator.append_lens_block(**l_1x2mm   )
collimating_transfocator.append_lens_block(**l_1x1mm   )
collimating_transfocator.append_lens_block(**l_1x500um )
collimating_transfocator.append_lens_block(**l_1x500um )
collimating_transfocator.append_lens_block(**l_1x300um )
collimating_transfocator.append_lens_block(**l_1x200um )
collimating_transfocator.append_lens_block(**l_2x200um )



print(collimating_transfocator.info())
print(collimating_transfocator.get_focal_distance())
print(collimating_transfocator.to_dictionary())

Fs = collimating_transfocator.calculate_focal_length_for_all_lens_block_status()
print(Fs.min(), Fs.max())

isorted = numpy.argsort(Fs)
x = numpy.arange(Fs.size)
plot(x, Fs[isorted], ylog=True)



photon_energy=10000.0
wanted_focal_distances = numpy.concatenate((numpy.linspace(0.1, 1000, 100), [numpy.inf]))
obtained_focal_distances = []
block_status = []

for focal_distance in wanted_focal_distances:
    collimating_transfocator.guess_and_set_lens_block_status(photon_energy=photon_energy, p=focal_distance)
    # collimating_transfocator.get_mask_from_index(idx)
    f_with_mask = collimating_transfocator.get_focal_distance(photon_energy=photon_energy)
    obtained_focal_distances.append(f_with_mask)
    mask = collimating_transfocator.get_lens_block_status()
    block_status.append(mask)
    print(">>>>>>--->", mask, f_with_mask)

plot(numpy.array(wanted_focal_distances), numpy.array(wanted_focal_distances),
     numpy.array(wanted_focal_distances), numpy.array(obtained_focal_distances),
     legend=["wanted","obtained"], xtitle="f wanted [m]", ytitle="f [m]")
#
for i in range(wanted_focal_distances.size):
    print(">>>", wanted_focal_distances[i], obtained_focal_distances[i], block_status[i])


# collimating_transfocator.set_lens_block_status([1,1,1,1,1,1,1,1,1])
print(collimating_transfocator.info_configuration(photon_energy=10000))

