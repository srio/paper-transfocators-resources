#!/usr/bin/python
# coding: utf-8
###################################################################################
# ID18 beamline - E=7 keV
# Authors/Contributors: Rafael Celestre
# Rafael.Celestre@esrf.eu
# creation: 01/04/2021
# last update: 14/04/2021 (v0.1)
###################################################################################

import array
import datetime
import logging.handlers
import numpy as np
import os
import sys
import time

sys.path.insert(0, './srw_python')
from srwlib import *
from uti_plot import *

from skimage.restoration import unwrap_phase

def get_radii(_wfr, stvt_x=50, stvt_y=50, silent=True):
    k = 2 * np.pi / srwl_uti_ph_en_conv(_wfr.mesh.eStart, _in_u='eV', _out_u='m')
    arP1 = array('d', [0] * _wfr.mesh.nx * _wfr.mesh.ny)
    srwl.CalcIntFromElecField(arP1, _wfr, 0, 4, 3, _wfr.mesh.eStart, 0, 0)

    wp_phase = np.reshape(arP1, (_wfr.mesh.ny, _wfr.mesh.nx))
    wp_phase_x = wp_phase[int(_wfr.mesh.ny / 2), int(_wfr.mesh.nx / 2) - stvt_x:int(_wfr.mesh.nx / 2) + stvt_x]
    wp_phase_y = wp_phase[int(_wfr.mesh.ny / 2) - stvt_y:int(_wfr.mesh.ny / 2) + stvt_y, int(_wfr.mesh.nx / 2)]

    uwp_phase_x = unwrap_phase(wp_phase_x)
    uwp_phase_y = unwrap_phase(wp_phase_y)

    dx = (_wfr.mesh.xFin - _wfr.mesh.xStart) / _wfr.mesh.nx
    dy = (_wfr.mesh.yFin - _wfr.mesh.yStart) / _wfr.mesh.ny

    nx = wp_phase_x.shape[0]
    ny = wp_phase_y.shape[0]

    xStart = - (dx * (nx - 1)) / 2.0
    xFin = xStart + dx * (nx - 1)

    yStart = - (dy * (ny - 1)) / 2.0
    yFin = yStart + dy * (ny - 1)

    x = np.linspace(xStart, xFin, nx)
    y = np.linspace(yStart, yFin, ny)

    px = np.polynomial.polynomial.polyfit(x, uwp_phase_x, 5)
    Rx = k / (2 * px[2])

    py = np.polynomial.polynomial.polyfit(y, uwp_phase_y, 5)
    Ry = k / (2 * py[2])

    if silent is False:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(3, 2)
        axs[0, 0].set_title("wrapped phase")
        im = axs[0, 0].imshow(wp_phase, extent=[_wfr.mesh.xStart * 1e6, _wfr.mesh.xFin * 1e6, _wfr.mesh.yStart * 1e6,
                                                _wfr.mesh.yFin * 1e6], cmap=plt.cm.binary_r)
        plt.colorbar(im, ax=axs[0, 0])

        axs[0, 1].set_title("unwrapped phase")
        im = axs[0, 1].imshow(unwrap_phase(wp_phase), extent=[_wfr.mesh.xStart * 1e6, _wfr.mesh.xFin * 1e6, _wfr.mesh.yStart * 1e6,
                                                 _wfr.mesh.yFin * 1e6], cmap=plt.cm.jet)
        plt.colorbar(im, ax=axs[0, 1])

        axs[1, 0].set_title("wrapped phase - fit")
        im = axs[1, 0].plot(x * 1e6, wp_phase_x, label='h')
        im = axs[1, 0].plot(y * 1e6, wp_phase_y, label='v')
        axs[1, 0].legend(loc=1)

        axs[1, 1].set_title("unwrapped phase")
        im = axs[1, 1].plot(x * 1e6, uwp_phase_x, label='h')
        im = axs[1, 1].plot(y * 1e6, uwp_phase_y, label='v')
        axs[1, 1].legend(loc=1)

        # Reconstructed phase
        ph_x = px[0] + px[1] * x + px[2] * x ** 2
        ph_y = py[0] + py[1] * y + py[2] * y ** 2

        axs[2, 0].set_title("reconstructed phase")
        im = axs[2, 0].plot(x * 1e6, ph_x, label='h')
        im = axs[2, 0].plot(y * 1e6, ph_y, label='v')
        axs[2, 0].legend(loc=1)

        axs[2, 1].set_title("residues")
        im = axs[2, 1].plot(x * 1e6, uwp_phase_x - ph_x, label='h')
        im = axs[2, 1].plot(y * 1e6, uwp_phase_y - ph_y, label='v')
        axs[2, 1].legend(loc=1)

        fig.tight_layout()
        plt.show()

    return Rx, Ry


if __name__=='__main__':
    startTime = time.time()
    print('###################') if (srwl_uti_proc_is_master()) else 0
    print('# ID18 simulation #') if (srwl_uti_proc_is_master()) else 0
    print('###################') if (srwl_uti_proc_is_master()) else 0

    #############################################################################
    # Program variables
    beamE = 7
    wfr_resolution = (512, 512)  # nx, ny
    screen_range = (-1E-3, 1E-3, -1E-3, 1E-3)  # x_Start, x_Fin, y_Start, y_Fin
    defocus = 0  # detecor position
    sampling_factor = 0  # sampling factor for adjusting nx, ny (effective if > 0)
    MultiE = False  # partially coherent simulation
    calculation = 0  # radiation characteristic to calculate (for MultiE):
    #   0- Total Intensity (s0);
    #   1- Four Stokes components;
    #   2- Mutual Intensity Cut vs X;
    #   3- Mutual Intensity Cut vs Y;
    #   4- Mutual Intensity Cuts and Degree of Coherence vs X & Y;
    #  10- Flux
    #  20- Electric Field (sum of fields from all macro-electrons, assuming CSR)
    #  40- Total Intensity, Mutual Intensity Cuts and Degree of Coherence vs X & Y;
    save = False
    plots = True
    nMacroElec = 5000  # total number of macro-electrons
    directory = 'srw_results'
    prfx = 'id18_c01_'
    strDataFolderName = directory

    energy = str(beamE)
    energy = energy.replace('.', 'p')

    position = str(defocus * 1e3)
    position = position.replace('.', 'p')

    if calculation == 0:
        calc = 'intensity'
    elif calculation == 2:
        calc = 'MIC_x'
    elif calculation == 3:
        calc = 'MIC_y'

    strIntPropOutFileName = prfx + energy + 'keV' + '_intensity.dat'
    strPhPropOutFileName  = prfx + energy + 'keV' + '_intensity.dat'
    strIntPrtlChrnc       = prfx + energy + 'keV_' + str(int(nMacroElec/1000))+'k_ME_' + calc + '.dat'

    #############################################################################
    #############################################################################
    # Logging all logging.infos

    # Get time stamp
    start0 = time.time()
    dt = datetime.datetime.fromtimestamp(start0).strftime('%Y-%m-%d_%H:%M:%S')

    # Initializing logging
    log = logging.getLogger('')
    log.setLevel(logging.INFO)
    format = logging.Formatter('%(levelname)s: %(message)s')

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(format)
    log.addHandler(ch)
    if (srwl_uti_proc_is_master()):
        fh = logging.handlers.RotatingFileHandler(os.path.join(os.getcwd(), strDataFolderName) + '/' + dt + '_'
            + strIntPropOutFileName.replace('_intensity.dat', '.log'), maxBytes=(1048576 * 5), backupCount=7)
        fh.setFormatter(format)
        log.addHandler(fh)

    #############################################################################
    #############################################################################
    # Photon source

    wavelength = srwl_uti_ph_en_conv(beamE, _in_u='keV', _out_u='m')
    k = 2*np.pi/wavelength
    z = 36.0

    # ******************************** Undulator parameters (CPMU18)
    numPer = 138	    # Number of ID Periods
    undPer = 0.018		# Period Length [m]
    phB = 0	        	# Initial Phase of the Horizontal field component
    sB = 1		        # Symmetry of the Horizontal field component vs Longitudinal position
    xcID = 0 			# Transverse Coordinates of Undulator Center [m]
    ycID = 0
    zcID = 0
    n = 1
    # ******************************** Storage ring parameters
    eBeam = SRWLPartBeam()
    eBeam.Iavg = 0.2             # average Current [A]
    eBeam.partStatMom1.x = 0
    eBeam.partStatMom1.y = 0
    eBeam.partStatMom1.z = -0.5*undPer*(numPer + 4)    # initial Longitudinal Coordinate (set before the ID)
    eBeam.partStatMom1.xp = 0  					       # initial Relative Transverse Velocities
    eBeam.partStatMom1.yp = 0

    # e- beam paramters (RMS) EBS
    sigEperE = 0.00089  # relative RMS energy spread
    sigX = 2.973e-05  # horizontal RMS size of e-beam [m]
    sigXp = 5.29e-06  # horizontal RMS angular divergence [rad]
    sigY = 4.37e-06  # vertical RMS size of e-beam [m]
    sigYp = 1.89e-06  # vertical RMS angular divergence [rad]
    eBeam.partStatMom1.gamma = 6.00 / 0.51099890221e-03  # Relative Energy

    n = 1
    if (2 * (2 * n * wavelength * eBeam.partStatMom1.gamma ** 2 / undPer - 1)) <= 0:
        n=3
        if (2 * (2 * n * wavelength * eBeam.partStatMom1.gamma ** 2 / undPer - 1)) <= 0:
            n = 5
            if (2 * (2 * n * wavelength * eBeam.partStatMom1.gamma ** 2 / undPer - 1)) <= 0:
                n = 7

    K = np.sqrt(2 * (2 * n * wavelength * eBeam.partStatMom1.gamma ** 2 / undPer - 1))
    B = K / (undPer * 93.3728962)  # Peak Horizontal field [T] (undulator)

    # 2nd order stat. moments
    eBeam.arStatMom2[0] = sigX*sigX			 # <(x-<x>)^2>
    eBeam.arStatMom2[1] = 0					 # <(x-<x>)(x'-<x'>)>
    eBeam.arStatMom2[2] = sigXp*sigXp		 # <(x'-<x'>)^2>
    eBeam.arStatMom2[3] = sigY*sigY		     # <(y-<y>)^2>
    eBeam.arStatMom2[4] = 0					 # <(y-<y>)(y'-<y'>)>
    eBeam.arStatMom2[5] = sigYp*sigYp		 # <(y'-<y'>)^2>
    eBeam.arStatMom2[10] = sigEperE*sigEperE # <(E-<E>)^2>/<E>^2

    # Electron trajectory
    eTraj = 0

    # Precision parameters
    arPrecSR = [0]*7
    arPrecSR[0] = 1		# SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    arPrecSR[1] = 0.01	# relative precision
    arPrecSR[2] = 0		# longitudinal position to start integration (effective if < zEndInteg)
    arPrecSR[3] = 0		# longitudinal position to finish integration (effective if > zStartInteg)
    arPrecSR[4] = 20000	# Number of points for trajectory calculation
    arPrecSR[5] = 1		# Use "terminating terms"  or not (1 or 0 respectively)
    arPrecSR[6] = sampling_factor # sampling factor for adjusting nx, ny (effective if > 0)
    sampFactNxNyForProp = arPrecSR[6] # sampling factor for adjusting nx, ny (effective if > 0)

    und = SRWLMagFldU([SRWLMagFldH(n, 'v', B, phB, sB, 1)], undPer, numPer)

    magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID]))

    # Monochromatic wavefront
    wfr = SRWLWfr()
    wfr.allocate(1, wfr_resolution[0], wfr_resolution[1])  # Photon Energy, Horizontal and Vertical Positions
    wfr.mesh.zStart = z
    wfr.mesh.eStart = beamE * 1E3
    wfr.mesh.xStart = screen_range[0]
    wfr.mesh.xFin = screen_range[1]
    wfr.mesh.yStart = screen_range[2]
    wfr.mesh.yFin = screen_range[3]
    wfr.partBeam = eBeam
    meshPartCoh = deepcopy(wfr.mesh)

    #############################################################################
    #############################################################################
    # Wavefront generation
    if srwl_uti_proc_is_master():
        # ******************************** Calculating Initial Wavefront and extracting Intensity:
        logging.info('- Performing Initial Electric Field calculation ... ')
        srwl.CalcElecFieldSR(wfr, eTraj, magFldCnt, arPrecSR)

        Rx, Ry = get_radii(wfr, stvt_x=50, stvt_y=50, silent=True)

        logging.info('Initial wavefront:')
        logging.info('Nx = %d, Ny = %d' % (wfr.mesh.nx, wfr.mesh.ny))
        logging.info('dx = %.4f um, dy = %.4f um' % ((wfr.mesh.xFin - wfr.mesh.xStart) * 1E6 / wfr.mesh.nx, (wfr.mesh.yFin - wfr.mesh.yStart) * 1E6 / wfr.mesh.ny))
        logging.info('range x = %.4f mm, range y = %.4f mm' % ((wfr.mesh.xFin - wfr.mesh.xStart) * 1E3, (wfr.mesh.yFin - wfr.mesh.yStart) * 1E3))
        logging.info('- Wavefront curvature:')
        logging.info('SRW native calculation: Rx = %.6f, Ry = %.6f' % (wfr.Rx, wfr.Ry))
        logging.info('Phase fit: Rx = %.6f, Ry = %.6f' % (Rx, Ry))
        logging.info('dRx = %.3f %%, dRy = %.3f %%' % ((Rx - wfr.Rx) * 100 / Rx, (Ry - wfr.Ry) * 100 / Ry))

    #############################################################################
    #############################################################################
    # Beamline assembly
    if srwl_uti_proc_is_master():
        logging.info('Setting up beamline')

        # ============= Coherence slit
        oeCoehSlit = SRWLOptA(_shape='r', _ap_or_ob='a', _Dx=40e-6, _Dy=206.8e-6)
        # ============= drift space
        oeD1 = SRWLOptD(65-z)
        # ============= first lens
        oeL1 = SRWLOptL(_Fx=41.10, _Fy=49.14)
        # ============= drift space
        oeD2 = SRWLOptD(170-65)
        # ============= second lens
        oeL2 = SRWLOptL(_Fx=26.21, _Fy=41.39)
        # ============= drift space
        oeD3= SRWLOptD(200-170)

        # ============= Wavefront Propagation Parameters =======================#
        #               [ 0] [1] [2]  [3]  [4]  [5]  [6]  [7]   [8]  [9] [10] [11]
        # OASYS parameters
        # ppCoehSlit =    [0,   0,  1,   0,   0,   1,   1,   1,    1,   0,   0,   0]
        # ppD1       =    [0,   0,  1,   3,   0,   1,   1,   1,    1,   0,   0,   0]
        # ppL1       =    [0,   0,  1,   0,   0,   1,   1,   1,    1,   0,   0,   0]
        # ppD2       =    [0,   0,  1,   4,   0,   3,   1,   3,    1,   0,   0,   0]
        # ppL2       =    [0,   0,  1,   0,   0,   1,   1,   1,    1,   0,   0,   0]
        # ppD3       =    [0,   0,  1,   4,   0,   2,   1,   2,    1,   0,   0,   0]

        # ============= Wavefront Propagation Parameters =======================#
        #
        ppCoehSlit =    [0,   0,  1,   0,   0,   1,   1,   1,    1,   0,   0,   0]
        ppD1       =    [0,   0,  1,   1,   0,   1,   4,   1,    4,   0,   0,   0]
        ppL1       =    [0,   0,  1,   0,   0,   1,   1,   1,    1,   0,   0,   0]
        ppD2       =    [0,   0,  1,   1,   0,   2,  .5,   2,   .5,   0,   0,   0]
        ppL2       =    [0,   0,  1,   0,   0,   1,   1,   1,    1,   0,   0,   0]
        ppD3       =    [0,   0,  1,   4,   0,   3,   1,   3,    1,   0,   0,   0]
        ppFinal    =    [0,   0,  1,   0,   0,  .1,   1,   0.1,  1,   0,   0,   0]


        '''
        [ 3]: Type of Free-Space Propagator:
               0- Standard Fresnel
               1- Fresnel with analytical treatment of the quadratic (leading) phase terms
               2- Similar to 1, yet with different processing near a waist
               3- For propagation from a waist over a ~large distance
               4- For propagation over some distance to a waist
        [ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
        [ 6]: Horizontal Resolution modification factor at Resizing
        [ 7]: Vertical Range modification factor at Resizing
        [ 8]: Vertical Resolution modification factor at Resizing
        '''

        srw_oe_array = []
        srw_pp_array = []

        # ============= Coherence slit
        srw_oe_array.append(oeCoehSlit)
        srw_pp_array.append(ppCoehSlit)
        # ============= drift space
        srw_oe_array.append(oeD1)
        srw_pp_array.append(ppD1)
        # # ============= L1
        srw_oe_array.append(oeL1)
        srw_pp_array.append(ppL1)
        # ============= drift space
        srw_oe_array.append(oeD2)
        srw_pp_array.append(ppD2)
        # # ============= L2
        srw_oe_array.append(oeL2)
        srw_pp_array.append(ppL2)
        # # ============= drift space
        srw_oe_array.append(oeD3)
        srw_pp_array.append(ppD3)
        # # ============= Recrop and adjust
        srw_pp_array.append(ppFinal)

        optBL = SRWLOptC(srw_oe_array, srw_pp_array)

    #############################################################################
    #############################################################################
    if (srwl_uti_proc_is_master()):
        # Electric field propagation
        logging.info('- Simulating Electric Field Wavefront Propagation ... ')
        srwl.PropagElecField(wfr, optBL)

        # TODO: compare with python FFT
        # srwl.SetRepresElecField(wfr, 'a')
        # wfr.unitElFldAng = 1

        logging.info('Propagated wavefront:')
        logging.info('Nx = %d, Ny = %d' % (wfr.mesh.nx, wfr.mesh.ny))
        logging.info('dx = %.4f um, dy = %.4f um' % ((wfr.mesh.xFin-wfr.mesh.xStart)*1E6/wfr.mesh.nx,
                                                     (wfr.mesh.yFin-wfr.mesh.yStart)*1E6/wfr.mesh.ny))
        logging.info('range x = %.4f um, range y = %.4f um' % ((wfr.mesh.xFin-wfr.mesh.xStart)*1E6,
                                                               (wfr.mesh.yFin-wfr.mesh.yStart)*1E6))
        logging.info('Rx = %.10f, Ry = %.10f' % (wfr.Rx, wfr.Ry))

        if save is True or plots is True:
            arI1 = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" 2D array to take intensity data
            srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
            arP1 = array('d', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" array to take 2D phase data (note it should be 'd')
            srwl.CalcIntFromElecField(arP1, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0)

        if save:
            srwl_uti_save_intens_ascii(arI1, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName,  strIntPropOutFileName), 0)
            srwl_uti_save_intens_ascii(arP1, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName, strPhPropOutFileName),0)

        logging.info('>> single electron calculations: done')


    if MultiE is True:
        logging.info('- Simulating Partially-Coherent Wavefront Propagation... ') if(srwl_uti_proc_is_master()) else 0
        nMacroElecAvgPerProc = 10   # number of macro-electrons / wavefront to average on worker processes
        nMacroElecSavePer = 100     # intermediate data saving periodicity (in macro-electrons)
        srCalcMeth = 1              # SR calculation method
        srCalcPrec = 0.01           # SR calculation rel. accuracy
        radStokesProp = srwl_wfr_emit_prop_multi_e(eBeam, magFldCnt, meshPartCoh, srCalcMeth, srCalcPrec, nMacroElec,
                                                   nMacroElecAvgPerProc, nMacroElecSavePer, os.path.join(os.getcwd(),
                                                   strDataFolderName, strIntPrtlChrnc), sampFactNxNyForProp, optBL, _char=calculation)
        logging.info('>> multi electron electron calculations: done') if(srwl_uti_proc_is_master()) else 0

    deltaT = time.time() - startTime
    hours, minutes = divmod(deltaT, 3600)
    minutes, seconds = divmod(minutes, 60)
    logging.info(">>>> Elapsed time: " + str(int(hours)) + "h " + str(int(minutes)) + "min " + str(seconds) + "s ") if (srwl_uti_proc_is_master()) else 0


    if plots is True:
        # ********************************Electrical field intensity and phase after propagation
        plotMesh1x = [1E6 * wfr.mesh.xStart, 1E6 * wfr.mesh.xFin, wfr.mesh.nx]
        plotMesh1y = [1E6 * wfr.mesh.yStart, 1E6 * wfr.mesh.yFin, wfr.mesh.ny]
        uti_plot2d(arI1, plotMesh1x, plotMesh1y,['Horizontal Position [um]', 'Vertical Position [um]', 'Intensity After Propagation'])
        uti_plot2d(arP1, plotMesh1x, plotMesh1y,['Horizontal Position [um]', 'Vertical Position [um]', 'Phase After Propagation'])

        uti_plot_show()