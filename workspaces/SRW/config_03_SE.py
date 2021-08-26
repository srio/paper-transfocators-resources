#!/usr/bin/python
# coding: utf-8
###################################################################################
# ID18 beamline - E=7 keV
# Authors/Contributors: Rafael Celestre
# Rafael.Celestre@esrf.eu
# creation: 01/04/2021
# last update: 19/07/2021 (v0.2)
###################################################################################

import sys
import os
import time

# try:
#     from oasys_srw.srwlib import *
#     from oasys_srw.uti_plot import *
# except:

sys.path.insert(0, './srw_python')
from srwlib import *
from uti_plot import *

####################################################
# DEBUG TOOLS

if __name__ == '__main__':
    startTime = time.time()
    print('##################################') if (srwl_uti_proc_is_master()) else 0
    print('# transfocators paper: config 03 #') if (srwl_uti_proc_is_master()) else 0
    print('##################################') if (srwl_uti_proc_is_master()) else 0

    plots = False
    save = True
    MultiE = True  # partially coherent simulation
    beamE = 7.0
    calculation = 0
    nMacroElec = 100000

    strDataFolderName = 'me_tests'
    prfx = 'id18_c03_'
    energy = str(beamE)
    strIntPropOutFileName = prfx + energy + 'keV' + '_intensity.dat'
    strPhPropOutFileName  = prfx + energy + 'keV' + '_intensity.dat'
    strIntPrtlChrnc       = prfx + energy + 'keV_' + str(int(nMacroElec/1000)).replace('.', 'p')+'k_ME_intensity.dat'

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
    magnetic_fields.append(SRWLMagFldH(1, 'v', _B=1.1013641049909946, _ph=0.0, _s=-1, _a=1.0))
    magnetic_structure = SRWLMagFldU(_arHarm=magnetic_fields, _per=0.018, _nPer=138.0)
    magnetic_field_container = SRWLMagFldC(_arMagFld=[magnetic_structure],
                                           _arXc=array('d', [0.0]),
                                           _arYc=array('d', [0.0]),
                                           _arZc=array('d', [0.0]))

    mesh = SRWLRadMesh(_eStart=beamE*1e3,
                       _eFin  =beamE*1e3,
                       _ne    =1,
                       _xStart=-0.001,
                       _xFin  =0.001,
                       _nx    =100,
                       _yStart=-0.001,
                       _yFin  =0.001,
                       _ny    =100,
                       _zStart=36.0)

    # stk = SRWLStokes()
    # stk.allocate(1,100,100)
    # stk.mesh = mesh
    #
    wfr = SRWLWfr()
    wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
    wfr.mesh = mesh
    wfr.partBeam = part_beam
    wfr.unitElFld = 1

    initial_mesh = deepcopy(wfr.mesh)

    meshPartCoh = deepcopy(wfr.mesh)

    if srwl_uti_proc_is_master():
        # ******************************** Calculating Initial Wavefront and extracting Intensity:
        print('- Performing Initial Electric Field calculation ... ')
        srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [1, 0.01, 0.0, 0.0, 50000, 1, 0.2])

        print('Initial wavefront:')
        print('Nx = %d, Ny = %d' % (wfr.mesh.nx, wfr.mesh.ny))
        print('dx = %.4f um, dy = %.4f um' % ((wfr.mesh.xFin - wfr.mesh.xStart) * 1E6 / wfr.mesh.nx, (wfr.mesh.yFin - wfr.mesh.yStart) * 1E6 / wfr.mesh.ny))
        print('range x = %.4f mm, range y = %.4f mm' % ((wfr.mesh.xFin - wfr.mesh.xStart) * 1E3, (wfr.mesh.yFin - wfr.mesh.yStart) * 1E3))
        print('- Wavefront curvature:')
        print('SRW native calculation: Rx = %.6f, Ry = %.6f' % (wfr.Rx, wfr.Ry))

    # if plots:
    #     mesh0 = deepcopy(wfr.mesh)
    #     arI = array('f', [0]*mesh0.nx*mesh0.ny)
    #     srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
    #     arIx = array('f', [0]*mesh0.nx)
    #     srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
    #     arIy = array('f', [0]*mesh0.ny)
    #     srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)
    #
    #     arP = array('d', [0]*mesh0.nx*mesh0.ny)
    #     srwl.CalcIntFromElecField(arP, wfr, 0, 4, 3, mesh0.eStart, 0, 0)
    #     arPx = array('d', [0]*mesh0.nx)
    #     srwl.CalcIntFromElecField(arPx, wfr, 0, 4, 1, mesh0.eStart, 0, 0)
    #     arPy = array('d', [0]*mesh0.ny)
    #     srwl.CalcIntFromElecField(arPy, wfr, 0, 4, 2, mesh0.eStart, 0, 0)
    #
    #     #save ascii file with intensity
    #     #srwl_uti_save_intens_ascii(arI, mesh0, <file_path>)
    #     plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
    #     plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
    #     uti_plot2d1d(arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])
    #     uti_plot2d1d(arP, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Phase Before Propagation'])

    ####################################################
    # BEAMLINE
    srw_oe_array = []
    srw_pp_array = []

    # ============= Coherence slit
    oe_0=SRWLOptA(_shape='r',
                   _ap_or_ob='a',
                   _Dx=85.1e-6,
                   _Dy=506.7e-6,
                   _x=0.0,
                   _y=0.0)

    pp_oe_0 = [0,0,1.0,1,0,1.0,5.0,1.0,5.0,0,0.0,0.0]

    srw_oe_array.append(oe_0)
    srw_pp_array.append(pp_oe_0)

    # ============= drift space
    drift_before_oe_1 = SRWLOptD(29.)
    pp_drift_before_oe_1 = [0, 0, 1, 1, 0, 1., 1., 1.4142, 1.4142, 0., 0., 0.]

    srw_oe_array.append(drift_before_oe_1)
    srw_pp_array.append(pp_drift_before_oe_1)

    # ============= transfocator 1
    Rh = 46.1 * 2 * 6.960764234342776e-06
    oe_2=srwl_opt_setup_CRL(_foc_plane=1,
                    _delta=6.960764234342776e-06,
                    _atten_len=0.0033406302038437438,
                    _shape=1,
                    _apert_h=0.0011,
                    _apert_v=0.0011,
                    _r_min=Rh,
                    _n=1,
                    _wall_thick=5e-05,
                    _xc=0.0,
                    _yc=0.0,
                    _void_cen_rad=None,
                    _e_start=6999.997980650263,
                    _e_fin=6999.997980650263,
                    _nx=1001,
                    _ny=1001)

    pp_oe_2 = [0, 0, 1, 1, 0, 1., 2., 1., 1., 0., 0., 0.]
    srw_oe_array.append(oe_2)
    srw_pp_array.append(pp_oe_2)

    Rv = 85.2 * 2 * 6.960764234342776e-06
    oe_3=srwl_opt_setup_CRL(_foc_plane=2,
                    _delta=6.960764234342776e-06,
                    _atten_len=0.0033406302038437438,
                    _shape=1,
                    _apert_h=0.0011,
                    _apert_v=0.0011,
                    _r_min=Rv,
                    _n=1,
                    _wall_thick=5e-05,
                    _xc=0.0,
                    _yc=0.0,
                    _void_cen_rad=None,
                    _e_start=6999.997980650263,
                    _e_fin=6999.997980650263,
                    _nx=1001,
                    _ny=1001)

    pp_oe_3 = [0, 0, 1, 1, 0, 1., 1., 1., 1., 0., 0., 0.]

    srw_oe_array.append(oe_3)
    srw_pp_array.append(pp_oe_3)


    oe_4=SRWLOptA(_shape='r',
                   _ap_or_ob='a',
                   _Dx=0.001,
                   _Dy=0.001,
                   _x=0.0,
                   _y=0.0)

    pp_oe_4 = [0, 0, 1, 0, 0, 1., 1., 1., 1., 0., 0., 0.]

    srw_oe_array.append(oe_4)
    srw_pp_array.append(pp_oe_4)

    # ============= drift space
    drift_before_oe_5 = SRWLOptD(105.0)
    # pp_drift_before_oe_5  = [0, 0, 1, 1, 0, 1.5, 1., 1., 2., 0., 0., 0.]    # zero padding to reduce aliasing
    pp_drift_before_oe_5  = [0, 0, 1, 1, 0, 1.4142, 1.4142, 1.4142, 1.4142, 0., 0., 0.]    # zero padding to reduce aliasing

    srw_oe_array.append(drift_before_oe_5)
    srw_pp_array.append(pp_drift_before_oe_5)

    # ============= transfocator 2
    Rh = 31.8 * 2 * 6.960764234342776e-06
    oe_6=srwl_opt_setup_CRL(_foc_plane=1,
                    _delta=6.960764234342776e-06,
                    _atten_len=0.0033406302038437438,
                    _shape=1,
                    _apert_h=0.0011,
                    _apert_v=0.0011,
                    _r_min=Rh,
                    _n=1,
                    _wall_thick=5e-05,
                    _xc=0.0,
                    _yc=0.0,
                    _void_cen_rad=None,
                    _e_start=6999.997980650263,
                    _e_fin=6999.997980650263,
                    _nx=1001,
                    _ny=1001)

    # pp_oe_6 = [0, 0, 1, 1, 0, 3/4, 4/3, 1/5, 5., 0., 0., 0.]
    pp_oe_6 = [0, 0, 1, 1, 0, 3/4, 4/3, 1, 1., 0., 0., 0.]

    srw_oe_array.append(oe_6)
    srw_pp_array.append(pp_oe_6)

    Rv = 27.8* 2 * 6.960764234342776e-06
    oe_7=srwl_opt_setup_CRL(_foc_plane=2,
                    _delta=6.960764234342776e-06,
                    _atten_len=0.0033406302038437438,
                    _shape=1,
                    _apert_h=0.0011,
                    _apert_v=0.0011,
                    _r_min=Rv,
                    _n=1,
                    _wall_thick=5e-05,
                    _xc=0.0,
                    _yc=0.0,
                    _void_cen_rad=None,
                    _e_start=6999.997980650263,
                    _e_fin=6999.997980650263,
                    _nx=1001,
                    _ny=1001)

    pp_oe_7 = [0, 0, 1, 1, 0, 1., 1., 1., 1., 0., 0., 0.]

    srw_oe_array.append(oe_7)
    srw_pp_array.append(pp_oe_7)
    #
    oe_8=SRWLOptA(_shape='r',
                   _ap_or_ob='a',
                   _Dx=0.001,
                   _Dy=0.001,
                   _x=0.0,
                   _y=0.0)

    pp_oe_8 = [0, 0, 1, 1, 0, 1., 1., 1., 1., 0., 0., 0.]

    srw_oe_array.append(oe_8)
    srw_pp_array.append(pp_oe_8)

    # ============= drift space to image plane
    drift_before_oe_9 = SRWLOptD(30)
    pp_drift_before_oe_9 = [0, 0, 1, 1, 0, 1., 1., 1.5, 1., 0., 0., 0.]

    #
    srw_oe_array.append(drift_before_oe_9)
    srw_pp_array.append(pp_drift_before_oe_9)

    # # # ============= reshape and resize
    pp_final = [0, 0, 1, 1, 0, 1/4, 1., 1/8, 1.5, 0., 0., 0.]
    srw_pp_array.append(pp_final)

    optBL = SRWLOptC(srw_oe_array, srw_pp_array)

    ####################################################
    # PROPAGATION
    if srwl_uti_proc_is_master():
        print('- Simulating Electric Field Wavefront Propagation ... ')
        srwl.PropagElecField(wfr, optBL)

        print('Initial wavefront:')
        print('Nx = %d, Ny = %d' % (wfr.mesh.nx, wfr.mesh.ny))
        print('dx = %.4f um, dy = %.4f um' % (
        (wfr.mesh.xFin - wfr.mesh.xStart) * 1E6 / wfr.mesh.nx, (wfr.mesh.yFin - wfr.mesh.yStart) * 1E6 / wfr.mesh.ny))
        print('range x = %.4f mm, range y = %.4f mm' % (
        (wfr.mesh.xFin - wfr.mesh.xStart) * 1E3, (wfr.mesh.yFin - wfr.mesh.yStart) * 1E3))
        print('- Wavefront curvature:')
        print('SRW native calculation: Rx = %.6f, Ry = %.6f' % (wfr.Rx, wfr.Ry))

        if save is True or plots is True:
            mesh = deepcopy(wfr.mesh)
            arI = array('f', [0]*mesh.nx*mesh.ny)
            srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh.eStart, 0, 0)
            arP = array('d', [0]*mesh.nx*mesh.ny)
            srwl.CalcIntFromElecField(arP, wfr, 0, 4, 3, mesh.eStart, 0, 0)
            if save:
                srwl_uti_save_intens_ascii(arI, wfr.mesh,
                                           os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName), 0)
                # srwl_uti_save_intens_ascii(arP, wfr.mesh,
                #                            os.path.join(os.getcwd(), strDataFolderName, strPhPropOutFileName), 0)
                print('>>> saved files')
        if plots:

            arIx = array('f', [0]*mesh.nx)
            srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh.eStart, 0, 0)
            arIy = array('f', [0]*mesh.ny)
            srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh.eStart, 0, 0)

            arPx = array('d', [0]*mesh.nx)
            srwl.CalcIntFromElecField(arPx, wfr, 0, 4, 1, mesh.eStart, 0, 0)
            arPy = array('d', [0]*mesh.ny)
            srwl.CalcIntFromElecField(arPy, wfr, 0, 4, 2, mesh.eStart, 0, 0)

            plotMeshx = [1000*mesh.xStart, 1000*mesh.xFin, mesh.nx]
            plotMeshy = [1000*mesh.yStart, 1000*mesh.yFin, mesh.ny]
            uti_plot2d1d(arI, plotMeshx, plotMeshy, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity After Propagation'])
            uti_plot2d1d(arP, plotMeshx, plotMeshy, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Phase After Propagation'])

            uti_plot_show()

    print('>> single electron calculations: done')

    if MultiE is True:
        print('- Simulating Partially-Coherent Wavefront Propagation... ') if(srwl_uti_proc_is_master()) else 0
        nMacroElecAvgPerProc = 10   # number of macro-electrons / wavefront to average on worker processes
        nMacroElecSavePer = 100     # intermediate data saving periodicity (in macro-electrons)
        srCalcMeth = 1              # SR calculation method
        srCalcPrec = 0.01           # SR calculation rel. accuracy
        radStokesProp = srwl_wfr_emit_prop_multi_e(part_beam,
                                                   magnetic_field_container,
                                                   meshPartCoh,
                                                   srCalcMeth,
                                                   srCalcPrec,
                                                   nMacroElec,
                                                   nMacroElecAvgPerProc,
                                                   nMacroElecSavePer,
                                                   os.path.join(os.getcwd(), strDataFolderName, strIntPrtlChrnc),
                                                   0.2,
                                                   optBL,
                                                   _char=calculation)
        print('>> multi electron electron calculations: done') if(srwl_uti_proc_is_master()) else 0


    deltaT = time.time() - startTime
    hours, minutes = divmod(deltaT, 3600)
    minutes, seconds = divmod(minutes, 60)
    print(">>>> Elapsed time: " + str(int(hours)) + "h " + str(int(minutes)) + "min " + str(seconds) + "s ") if (srwl_uti_proc_is_master()) else 0

