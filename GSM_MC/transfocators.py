import sys

from sr.crl import LensBlock, Transfocator

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
#        l_1x300um,
        l_1x200um,
        l_2x200um,
        l_4x200um,
#        l_8x200um,
    )
)

collimating_transfocator_h = Transfocator(
    (
#        l_1x10mm,
        l_1x5mm,
        l_1x2mm,
        l_1x1mm,
        l_1x500um,
        l_1x200um,
        l_2x200um,
    )
)


FO1 = Transfocator(
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

FO1_H = Transfocator(
    (
        l_1x5mm,
        l_1x2mm,
        l_1x1mm,
        l_1x500um,
#        l_1x300um,
#        l_1x200um,
#        l_2x200um,
    )
)


FO2 = Transfocator(
    (
#        l_1x5mm,
        l_1x2mm,
        l_1x1mm,
        l_1x500um,
 #       l_1x300um,
        l_1x200um,
        l_2x200um,
        l_4x200um,
        l_8x200um,
        l_16x200um,
        l_32x200um,
#        l_32x200um, # need for 35 keV in EH1
    )
)

FO2_H = Transfocator(
    (
        l_1x5mm,
        l_1x1mm,
        l_1x500um,
        l_1x300um,
#        l_1x200um,
#        l_2x200um,
    )
)



TRANSFOCATORS_2D = {
    65: collimating_transfocator,
    170: FO1,
    192: FO2,
}

TRANSFOCATORS_H = {
    65: collimating_transfocator_h,
    170: FO1_H,
    192: FO2_H,
}


for pos in TRANSFOCATORS_2D.keys():
    t = TRANSFOCATORS_2D[pos]+TRANSFOCATORS_H[pos]
    print(f"\n\n## TRANSFOCATOR @ {pos}")
    print(t)
