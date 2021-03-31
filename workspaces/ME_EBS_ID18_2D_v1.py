# run using OAR:
# oarsub -S ./oar.bash.sh
# watch oarstat -u

try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy

#if not srwl_uti_proc_is_master(): exit()

####################################################
# LIGHT SOURCE

part_beam = SRWLPartBeam()
part_beam.Iavg               = 0.2
part_beam.partStatMom1.x     = 0.0
part_beam.partStatMom1.y     = 0.0
part_beam.partStatMom1.z     = -1.3139999999999998
part_beam.partStatMom1.xp    = 0.0
part_beam.partStatMom1.yp    = 0.0
part_beam.partStatMom1.gamma = 11741.70710144324
part_beam.arStatMom2[0]      = 8.838729000000001e-10
part_beam.arStatMom2[1]      = 0.0
part_beam.arStatMom2[2]      = 1.9096899999999998e-11
part_beam.arStatMom2[3]      = 2.79841e-11
part_beam.arStatMom2[4]      = 0.0
part_beam.arStatMom2[5]      = 3.5720999999999997e-12
part_beam.arStatMom2[10]     = 7.920999999999999e-07

magnetic_fields = []
magnetic_fields.append(SRWLMagFldH(1, 'v', 
                                   _B=1.1013641049909946, 
                                   _ph=0.0, 
                                   _s=-1, 
                                   _a=1.0))
magnetic_structure = SRWLMagFldU(_arHarm=magnetic_fields, _per=0.018, _nPer=138.0)
magnetic_field_container = SRWLMagFldC(_arMagFld=[magnetic_structure], 
                                       _arXc=array('d', [0.0]), 
                                       _arYc=array('d', [0.0]), 
                                       _arZc=array('d', [0.0]))

mesh = SRWLRadMesh(_eStart=6999.997980650263,
                   _eFin  =6999.997980650263,
                   _ne    =1,
                   _xStart=-0.001,
                   _xFin  =0.001,
                   _nx    =100,
                   _yStart=-0.001,
                   _yFin  =0.001,
                   _ny    =100,
                   _zStart=36.0)

stk = SRWLStokes()
stk.allocate(1,100,100)
stk.mesh = mesh

wfr = SRWLWfr()
wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
wfr.mesh = mesh
wfr.partBeam = part_beam
wfr.unitElFld = 1

initial_mesh = deepcopy(wfr.mesh)

####################################################
# BEAMLINE

srw_oe_array = []
srw_pp_array = []

oe_0=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=9.4e-05,
               _Dy=0.000486,
               _x=0.0,
               _y=0.0)

pp_oe_0 = [0,0,1.0,0,0,1.0,20.0,1.0,20.0,0,0.0,0.0]

srw_oe_array.append(oe_0)
srw_pp_array.append(pp_oe_0)

drift_before_oe_1 = SRWLOptD(29.0)
pp_drift_before_oe_1 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_1)
srw_pp_array.append(pp_drift_before_oe_1)


oe_2=SRWLOptL(_Fx=41.04, _Fy=49.18, _x=0.0, _y=0.0)

pp_oe_2 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_2)
srw_pp_array.append(pp_oe_2)

oe_3=SRWLOptA(_shape='c',
               _ap_or_ob='a',
               _Dx=0.002,
               _Dy=0.002,
               _x=0.0,
               _y=0.0)

pp_oe_3 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_3)
srw_pp_array.append(pp_oe_3)

drift_before_oe_4 = SRWLOptD(105.0)
pp_drift_before_oe_4 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_4)
srw_pp_array.append(pp_drift_before_oe_4)


oe_5=SRWLOptL(_Fx=26.21, _Fy=41.28, _x=0.0, _y=0.0)

pp_oe_5 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_5)
srw_pp_array.append(pp_oe_5)

oe_6=SRWLOptA(_shape='c',
               _ap_or_ob='a',
               _Dx=2.0,
               _Dy=2.0,
               _x=0.0,
               _y=0.0)

pp_oe_6 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_6)
srw_pp_array.append(pp_oe_6)

drift_before_oe_7 = SRWLOptD(30.0)
pp_drift_before_oe_7 = [0,0,1.0,1,0,0.5,2.0,0.5,2.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_7)
srw_pp_array.append(pp_drift_before_oe_7)



####################################################
# PROPAGATION

optBL = SRWLOptC(srw_oe_array, srw_pp_array)



####################################################
# MULTI ELECTRON PROPAGATION

radStokesProp = srwl_wfr_emit_prop_multi_e(part_beam,
                                           magnetic_field_container,
                                           initial_mesh,
                                           1,
                                           0.01,
                                           500000,
                                           5,
                                           20,
                                           'output_srw_script_me.dat',
                                           1.0,
                                           optBL,
                                           _char=0)
