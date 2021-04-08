from oasys_srw.srwlib import *
import numpy as np

wfr = in_object_1._SRWData__srw_wavefront
k = 2*np.pi/srwl_uti_ph_en_conv(wfr.mesh.eStart, _in_u='eV', _out_u='m')

arP1 = array('d', [0] * wfr.mesh.nx * wfr.mesh.ny) 
srwl.CalcIntFromElecField(arP1, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0)
wp_phase = np.reshape(arP1, (wfr.mesh.ny, wfr.mesh.nx))

wp_phase_x = np.unwrap(wp_phase[int(wfr.mesh.ny/2),int(wfr.mesh.nx/2)-25:int(wfr.mesh.nx/2)+25])
wp_phase_y = np.unwrap(wp_phase[int(wfr.mesh.ny/2)-25:int(wfr.mesh.ny/2)+25,int(wfr.mesh.nx/2)])

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

p = np.polynomial.polynomial.polyfit(x,wp_phase_x,5)
Rx = k/(2*p[2])

p = np.polynomial.polynomial.polyfit(x,wp_phase_y,5)
Ry = k/(2*p[2])

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


