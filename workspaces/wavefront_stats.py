wfr = in_object_1._SRWData__srw_wavefront

print('Propagated wavefront:')
print('Nx = %d, Ny = %d' % (wfr.mesh.nx, wfr.mesh.ny))
print('dx = %.4f um, dy = %.4f um' % ((wfr.mesh.xFin-wfr.mesh.xStart)*1E6/wfr.mesh.nx,(wfr.mesh.yFin-wfr.mesh.yStart)*1E6/wfr.mesh.ny))
print('range x = %.4f um, range y = %.4f um' % ((wfr.mesh.xFin-wfr.mesh.xStart)*1E6,(wfr.mesh.yFin-wfr.mesh.yStart)*1E6))
print('Rx = %.10f, Ry = %.10f' % (wfr.Rx, wfr.Ry))