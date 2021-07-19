try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy

if not srwl_uti_proc_is_master(): exit()

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
srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [1,0.01,0.0,0.0,50000,1,0.1])

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

####################################################
# BEAMLINE

srw_oe_array = []
srw_pp_array = []

oe_0=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=4.03e-05,
               _Dy=0.000227,
               _x=0.0,
               _y=0.0)

pp_oe_0 = [0,0,1.0,1,0,1.0,10.0,1.0,10.0,0,0.0,0.0]

srw_oe_array.append(oe_0)
srw_pp_array.append(pp_oe_0)

drift_before_oe_1 = SRWLOptD(29.0)
pp_drift_before_oe_1 = [0,0,1.0,3,0,2.0,1.0,2.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_1)
srw_pp_array.append(pp_drift_before_oe_1)


oe_2=srwl_opt_setup_CRL(_foc_plane=1,
                _delta=6.960764234342776e-06,
                _atten_len=0.0033406302038437438,
                _shape=1,
                _apert_h=0.001,
                _apert_v=0.001,
                _r_min=0.0006418999999999999,
                _n=1,
                _wall_thick=2.9999999999999997e-05,
                _xc=0.0,
                _yc=0.0,
                _void_cen_rad=None,
                _e_start=6999.997980650263,
                _e_fin=6999.997980650263,
                _nx=5001,
                _ny=5001)

pp_oe_2 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_2)
srw_pp_array.append(pp_oe_2)

oe_3=srwl_opt_setup_CRL(_foc_plane=2,
                _delta=6.960764234342776e-06,
                _atten_len=0.0033406302038437438,
                _shape=1,
                _apert_h=0.001,
                _apert_v=0.001,
                _r_min=0.0002094,
                _n=1,
                _wall_thick=2.9999999999999997e-05,
                _xc=0.0,
                _yc=0.0,
                _void_cen_rad=None,
                _e_start=6999.997980650263,
                _e_fin=6999.997980650263,
                _nx=5001,
                _ny=5001)

pp_oe_3 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_3)
srw_pp_array.append(pp_oe_3)

oe_4=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=0.002,
               _Dy=0.002,
               _x=0.0,
               _y=0.0)

pp_oe_4 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_4)
srw_pp_array.append(pp_oe_4)

drift_before_oe_5 = SRWLOptD(105.0)
pp_drift_before_oe_5 = [0,0,1.0,1,0,2.0,1.0,2.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_5)
srw_pp_array.append(pp_drift_before_oe_5)


oe_6=srwl_opt_setup_CRL(_foc_plane=1,
                _delta=6.960764234342776e-06,
                _atten_len=0.0033406302038437438,
                _shape=1,
                _apert_h=0.001,
                _apert_v=0.001,
                _r_min=0.0003695,
                _n=1,
                _wall_thick=2.9999999999999997e-05,
                _xc=0.0,
                _yc=0.0,
                _void_cen_rad=None,
                _e_start=6999.997980650263,
                _e_fin=6999.997980650263,
                _nx=5001,
                _ny=5001)

pp_oe_6 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_6)
srw_pp_array.append(pp_oe_6)

oe_7=srwl_opt_setup_CRL(_foc_plane=2,
                _delta=6.960764234342776e-06,
                _atten_len=0.0033406302038437438,
                _shape=1,
                _apert_h=0.001,
                _apert_v=0.001,
                _r_min=0.0003096,
                _n=1,
                _wall_thick=2.9999999999999997e-05,
                _xc=0.0,
                _yc=0.0,
                _void_cen_rad=None,
                _e_start=6999.997980650263,
                _e_fin=6999.997980650263,
                _nx=5001,
                _ny=5001)

pp_oe_7 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_7)
srw_pp_array.append(pp_oe_7)

oe_8=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=0.002,
               _Dy=0.002,
               _x=0.0,
               _y=0.0)

pp_oe_8 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_8)
srw_pp_array.append(pp_oe_8)

drift_before_oe_9 = SRWLOptD(30.0)
pp_drift_before_oe_9 = [0,0,1.0,0,0,0.7,2.0,0.7,2.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_9)
srw_pp_array.append(pp_drift_before_oe_9)


oe_10=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=1.0,
               _Dy=1.0,
               _x=0.0,
               _y=0.0)

pp_oe_10 = [0,0,1.0,0,0,0.075,2.0,0.005,10.0,0,0.0,0.0]

srw_oe_array.append(oe_10)
srw_pp_array.append(pp_oe_10)


####################################################
# PROPAGATION

optBL = SRWLOptC(srw_oe_array, srw_pp_array)
srwl.PropagElecField(wfr, optBL)

mesh1 = deepcopy(wfr.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny)
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0)
arI1x = array('f', [0]*mesh1.nx)
srwl.CalcIntFromElecField(arI1x, wfr, 6, 0, 1, mesh1.eStart, 0, 0)
arI1y = array('f', [0]*mesh1.ny)
srwl.CalcIntFromElecField(arI1y, wfr, 6, 0, 2, mesh1.eStart, 0, 0)
#save ascii file with intensity
#srwl_uti_save_intens_ascii(arI1, mesh1, <file_path>)
plotMesh1x = [1000*mesh1.xStart, 1000*mesh1.xFin, mesh1.nx]
plotMesh1y = [1000*mesh1.yStart, 1000*mesh1.yFin, mesh1.ny]
uti_plot2d1d(arI1, plotMesh1x, plotMesh1y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity After Propagation'])
uti_plot_show()
