from srxraylib.plot.gol import plot, plot_image, plot_show
import numpy
import h5py
from tools import fwhm
import matplotlib.pyplot as plt


filename = "wofry1d_h_doubleslit.dat"
cf = numpy.loadtxt(filename)
APERTURES = cf[:,0]
CF = cf[:,1]
DoC = cf[:,2]


# plot(APERTURES*1e6,CF, xtitle="aperture(outer) [um]", ytitle="CF")
# plot(APERTURES*1e6,DoC, xtitle="aperture(outer) [um]", ytitle="DoC fitted")



f = h5py.File("wofry1d_h_doubleslit.h5", 'r')
SD = f["images/Intensity/image_data"][()].T
APERTURES = f["images/Intensity/axis_x"][()]
x = f["images/Intensity/axis_y"][()]
f.close()






ii = 9


###
plt.rcParams.update({'font.size': 28})
plot(x*1e6, SD[ii, :] / SD.max(), xtitle="$x$ [$\mu$m] ($z$=46 m)", ytitle="$\mathcal{I}$ [arbitrary units] ",
     xrange=[-1500,1500], title = "$s$ = %3.1f $\mu$m" % (APERTURES[ii]*1e6 - 2.5), show=0, figsize=(16, 12) )

filename = "doubleslit_profile.pdf"
plt.savefig(filename)
print("File written to disk: %s" % filename)
plot_show()

###

# SD[ii, :] = 0
plt.rcParams.update({'font.size': 20})
plot_image(SD / SD.max(), APERTURES*1e6 - 2.5, x*1e6, xtitle="$s$ [$\mu$m]",ytitle="$x$ [$\mu$m] ($z$=46 m)",
           title="", aspect='auto', xrange=[18,120], yrange=[-750,750], add_colorbar=False,figsize=(10,18), show=0)

filename = "doubleslit_scan.pdf"
plt.savefig(filename)
print("File written to disk: %s" % filename)
plot_show()


###
DoC_wofry1 = numpy.loadtxt("DoCprofile_wofry1.dat")
DoC_srwme = numpy.loadtxt("DoCprofile_srwme.dat")

APERTURES2 = []
DoC2 = []
i = 0
while i < APERTURES.size:
    APERTURES2.append(APERTURES[i])
    DoC2.append(DoC[i])
    i += 3
APERTURES2 = numpy.array(APERTURES2)
DoC2 = numpy.array(DoC2)

plt.rcParams.update({'font.size': 28})
plot(
    DoC_wofry1[:, 0], DoC_wofry1[:, 1],
    DoC_srwme[:,0], DoC_srwme[:,1],
    1e6 * APERTURES2 - 2.50, DoC2,
    title="", legend=[
        "WOFRY (CL: %3.1f $\mu$m) " %   (fwhm(DoC_srwme[:,0], DoC_srwme[:,1])),
        "SRW (CL: %3.1f $\mu$m) "   %   (fwhm(DoC_wofry1[:,0], DoC_wofry1[:,1])),
        "From $\mathcal{V}$ in double-slit patterns",
        ],
    marker=[None, None, '.'], linestyle=[None,None,""],
    xtitle="$x_2-x_1, s$ [$\mu$m]", ytitle="|DoC|", xrange=[0,150], yrange=[0, 1.1], show=0, figsize=(16,12))

filename = "doubleslit_DoC.pdf"
plt.savefig(filename)
print("File written to disk: %s" % filename)
plot_show()