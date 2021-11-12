try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy
from srxraylib.plot.gol import set_qt

if not srwl_uti_proc_is_master(): exit()

####################################################
# LIGHT SOURCE

part_beam = SRWLPartBeam()
part_beam.Iavg               = 0.2
part_beam.partStatMom1.x     = 0.0
part_beam.partStatMom1.y     = 0.0
part_beam.partStatMom1.z     = -0.8300000000000001
part_beam.partStatMom1.xp    = 0.0
part_beam.partStatMom1.yp    = 0.0
part_beam.partStatMom1.gamma = 11741.70710144324
part_beam.arStatMom2[0]      = 2.1904000000000003e-10
part_beam.arStatMom2[1]      = 0.0
part_beam.arStatMom2[2]      = 7.839999999999999e-12
part_beam.arStatMom2[3]      = 1.3690000000000002e-11
part_beam.arStatMom2[4]      = 0.0
part_beam.arStatMom2[5]      = 2.2500000000000003e-12
part_beam.arStatMom2[10]     = 1.9043999999999997e-06

magnetic_fields = []
magnetic_fields.append(SRWLMagFldH(1, 'v',
                                   _B=0.8032309541789979,
                                   _ph=0.0,
                                   _s=-1,
                                   _a=1.0))
magnetic_structure = SRWLMagFldU(_arHarm=magnetic_fields, _per=0.02, _nPer=75)
magnetic_field_container = SRWLMagFldC(_arMagFld=[magnetic_structure],
                                       _arXc=array('d', [0.0]),
                                       _arYc=array('d', [0.0]),
                                       _arZc=array('d', [0.0]))

mesh = SRWLRadMesh(_eStart=8043.9597627810945,
                   _eFin  =8043.9597627810945,
                   _ne    =1,
                   _xStart=-0.0005,
                   _xFin  =0.0005,
                   _nx    =100,
                   _yStart=-0.0005,
                   _yFin  =0.0005,
                   _ny    =100,
                   _zStart=10.0)

stk = SRWLStokes()
stk.allocate(1,100,100)
stk.mesh = mesh

wfr = SRWLWfr()
wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
wfr.mesh = mesh
wfr.partBeam = part_beam
wfr.unitElFld = 1

initial_mesh = deepcopy(wfr.mesh)
srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [1,0.01,0.0,0.0,50000,1,0.0])

mesh0 = deepcopy(wfr.mesh)
arI = array('f', [0]*mesh0.nx*mesh0.ny)
srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
arIx = array('f', [0]*mesh0.nx)
srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
arIy = array('f', [0]*mesh0.ny)
srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)
#save ascii file with intensity
#srwl_uti_save_intens_ascii(arI, mesh0, <file_path>)
plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
uti_plot2d1d (arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])
uti_plot_show()
