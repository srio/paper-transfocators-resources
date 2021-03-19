import numpy as np
import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.run as rr
import xrt.backends.raycing.materials as rm

mBeryllium = rm.Material('Be', rho=1.848, kind='lens')
mAluminum = rm.Material('Al', rho=2.7, kind='lens')

materials = dict()
materials["Be"]=mBeryllium
materials["Al"]=mAluminum

Lens = roe.DoubleParaboloidLens
from sr.crl import LensBlock, Transfocator

# lp = dict(aperture=1e-3)

lp = dict(thickness=1e-3)

l_1x10mm = LensBlock(1, radius=10000e-6, **lp)
l_1x5mm = LensBlock(1, radius=5000e-6, **lp)
l_1x2mm = LensBlock(1, radius=2000e-6, **lp)
l_1x1mm = LensBlock(1, radius=1000e-6, **lp)
l_1x500um = LensBlock(1, radius=500e-6, **lp)
l_1x500um = LensBlock(1, radius=500e-6, **lp)
l_1x300um = LensBlock(1, radius=300e-6, **lp)
l_1x200um = LensBlock(1, radius=200e-6, **lp)
l_2x200um = LensBlock(2, radius=200e-6, **lp)
l_4x200um = LensBlock(4, radius=200e-6, **lp)
l_8x200um = LensBlock(8, radius=200e-6, **lp)
l_16x200um = LensBlock(16, radius=200e-6, **lp)
l_32x200um = LensBlock(32, radius=200e-6, **lp)
l_64x200um = LensBlock(64, radius=200e-6, **lp)

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



class LensStack(Lens):
    def __init__(self,
            bl = None,
            name = "CRLstack",
            center = [0,30e3,0],
            material = "Be",
            radius = 1,
            thickness = 0.1,
            zmax = 1,
            pitch=np.pi/2,
            n=1):
        """
        radius and thickness in mm
        n is number of lens or distance,energy tuple
        """
        material = materials[material]
        Lens.__init__(self,
            bl=bl,
            name=name,
            material=material,
            # for xrt x^2/(4f), this is usually x^2/(2r), so f=r/2
            focus=radius/2,
            center=center,
            zmax=zmax,
            pitch=pitch,
            t=thickness,
            nCRL=n
        )
        if isinstance(n,(float,int)):
            _n = n
        else:
            _n = self.get_nCRL(*n)
            print("Lens stack would need %.1f lenses to focus at %.3f m at %.3f keV"%(_n,n[0]/1e3,n[1]/1e3))




class Transfocator:
    def __init__(self,
            sr_lensset,
            bl = None,
            name="Transfocator",
            center = [0,30e3,0],
            pitch=np.pi/2,
            ):
        """
        radius and thickness in mm
        n is number of lens or distance,energy tuple
        """
        _xrt_stacks = []
        # define XRT lens stacks
        for i,stack in enumerate(sr_lensset):
            radius = stack.radius*1e3 # sr uses m, xrt mm
            material = stack.material
            # for XRT thickness is at the center of the lens
            # in sr terminology this is web_thickness
            thickness = stack.web_thickness*1e3
            # for XRT double curved lenses the overall thickness is
            # 2*zmax+thickness
            zmax = (stack.thickness-stack.web_thickness)/2*1e3
            n = stack.n
            center = [center[0],center[1]+i*10*zmax,center[2]]

            xrt_stack = LensStack(
                bl=bl,
                name=f"{name}_crl_stack_{i}",
                material=material,
                radius=radius,
                center=center,
                zmax=zmax,
                pitch=pitch,
                thickness=thickness,
                n=n
            )
            _xrt_stacks.append(xrt_stack)
        self._xrt_stacks = _xrt_stacks

    def multiple_refract(self,beam=None):
        _1 = _2 = None # so that the return does no complaim if empty list
        for i,s in enumerate(self._xrt_stacks):
            print("transfocator propagation",i)
            beam,_1,_2=s.multiple_refract(beam)
        return beam,_1,_2

