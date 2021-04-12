from oasys_srw.srwlib import *
import numpy as np
from skimage.restoration import unwrap_phase
import matplotlib.pyplot as plt

# parameters

stvt_x = 8 # half window size in pixel for fit
stvt_y = 8 # half window size in pixel for fit
plots = False


wfr = in_object_1._SRWData__srw_wavefront
k = 2*np.pi/srwl_uti_ph_en_conv(wfr.mesh.eStart, _in_u='eV', _out_u='m')

# 2D phase extraction

arP1 = array('d', [0] * wfr.mesh.nx * wfr.mesh.ny) 
srwl.CalcIntFromElecField(arP1, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0)

wp_phase = np.reshape(arP1, (wfr.mesh.ny, wfr.mesh.nx))
wp_phase_x = wp_phase[int(wfr.mesh.ny/2),int(wfr.mesh.nx/2)-stvt_x:int(wfr.mesh.nx/2)+stvt_x]
wp_phase_y = wp_phase[int(wfr.mesh.ny/2)-stvt_y:int(wfr.mesh.ny/2)+stvt_y,int(wfr.mesh.nx/2)]

# Unwrapped phase
uwp_phase = unwrap_phase(wp_phase)
uwp_phase_x = unwrap_phase(wp_phase_x)
uwp_phase_y = unwrap_phase(wp_phase_y)

dx = (wfr.mesh.xFin - wfr.mesh.xStart) / wfr.mesh.nx
dy = (wfr.mesh.yFin - wfr.mesh.yStart) / wfr.mesh.ny

nx = wp_phase_x.shape[0]
ny = wp_phase_y.shape[0]

xStart = - (dx * (nx - 1)) / 2.0
xFin = xStart + dx * (nx - 1)

yStart = - (dy * (ny - 1)) / 2.0
yFin = yStart + dy * (ny - 1)

x = np.linspace(xStart, xFin, nx)
y = np.linspace(yStart, yFin, ny)

px = np.polynomial.polynomial.polyfit(x,uwp_phase_x,5)
Rx = k/(2*px[2])

py = np.polynomial.polynomial.polyfit(y,uwp_phase_y,5)
Ry = k/(2*py[2])

print()
print('- Propagated wavefront:')
print('Nx = %d, Ny = %d' % (wfr.mesh.nx, wfr.mesh.ny))
print('dx = %.4f um, dy = %.4f um' % ((wfr.mesh.xFin-wfr.mesh.xStart)*1E6/wfr.mesh.nx,(wfr.mesh.yFin-wfr.mesh.yStart)*1E6/wfr.mesh.ny))
print('range x = %.4f um, range y = %.4f um' % ((wfr.mesh.xFin-wfr.mesh.xStart)*1E6,(wfr.mesh.yFin-wfr.mesh.yStart)*1E6))
print()
print('- Wavefront curvature:')
print('SRW native: Rx = %.10f, Ry = %.10f' % (wfr.Rx, wfr.Ry))
print('Phase fit: Rx = %.10f, Ry = %.10f' % (Rx, Ry))
print('dRx = %.3f %%, dRy = %.3f %%' % ((Rx-wfr.Rx)*100/Rx, (Ry-wfr.Ry)*100/Ry))

if plots:
    # plots for visual inspection

    fig, axs = plt.subplots(3, 2)

    # arI1 = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny) 
    # srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
    # intensity = np.reshape(arP1, (wfr.mesh.ny, wfr.mesh.nx))

    axs[0, 0].set_title("wrapped phase")
    im = axs[0, 0].imshow(wp_phase, extent=[wfr.mesh.xStart*1e6, wfr.mesh.xFin*1e6, wfr.mesh.yStart*1e6, wfr.mesh.yFin*1e6],cmap=plt.cm.binary_r)
    plt.colorbar(im, ax=axs[0, 0])

    axs[0, 1].set_title("unwrapped phase")
    im = axs[0, 1].imshow(uwp_phase, extent=[wfr.mesh.xStart*1e6, wfr.mesh.xFin*1e6, wfr.mesh.yStart*1e6, wfr.mesh.yFin*1e6],cmap=plt.cm.jet)
    plt.colorbar(im, ax=axs[0, 1])

    axs[1, 0].set_title("wrapped phase - fit")
    im = axs[1, 0].plot(x*1e6, wp_phase_x, label='h')
    im = axs[1, 0].plot(y*1e6, wp_phase_y, label='v')
    axs[1, 0].legend(loc=1)

    axs[1, 1].set_title("unwrapped phase")
    im = axs[1, 1].plot(x*1e6, uwp_phase_x, label='h')
    im = axs[1, 1].plot(y*1e6, uwp_phase_y, label='v')
    axs[1, 1].legend(loc=1)

    # Reconstructed phase
    ph_x = px[0] + px[1]*x + px[2] * x**2
    ph_y = py[0] + py[1]*x + py[2] * x**2

    axs[2, 0].set_title("reconstructed phase")
    im = axs[2, 0].plot(x*1e6, ph_x, label='h')
    im = axs[2, 0].plot(y*1e6, ph_y, label='v')
    axs[2, 0].legend(loc=1)

    axs[2, 1].set_title("residues")
    im = axs[2, 1].plot(x*1e6, uwp_phase_x-ph_x, label='h')
    im = axs[2, 1].plot(y*1e6, uwp_phase_y-ph_y, label='v')
    axs[2, 1].legend(loc=1)

    fig.tight_layout()

    plt.subplot_tool()
    plt.show()


