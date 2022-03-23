#
# Import section
#
import numpy
from tools import plotCSD, plot_four_images_with_histograms, plot_show
from wofry1d_h_source import main as main_h
from wofry1d_v_source import main as main_v
import matplotlib.pylab as plt

#
# MAIN========================
#

from tools import plotCSD

tally_h = main_h(do_plot=0)


CSD_H, X_H, Y_H = plotCSD(
    tally_h, range_limits=[-100,100], compare_profiles=1, rotate_axes_flag=0, rotate_axes_normalization=2, normalize_to_DoC=0,
    direction='x',
    srw_file='/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/source_imaging/CSDx7.0keV_50k_ME_intensity.dat',
    do_plot=0)



tally_v = main_v(do_plot=0)

CSD_V, X_V, Y_V = plotCSD(
    tally_v, range_limits=[-100,100], compare_profiles=1, rotate_axes_flag=0, rotate_axes_normalization=2, normalize_to_DoC=0,
    direction='y',
    srw_file='/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/source_imaging/CSDy7.0keV_15k_ME_intensity.dat',
    do_plot=0)

CSD_H.append(CSD_V[0])
CSD_H.append(CSD_V[1])
X_H.append(X_V[0])
X_H.append(X_V[1])
Y_H.append(Y_V[0])
Y_H.append(Y_V[1])

plot_four_images_with_histograms(CSD_H, X_H, Y_H,
                                      title_list=["WOFRY1D H", "SRW-ME H", "WOFRY1D V", "SRW-ME V"],
                                      xtitle="$x$ [$\mu$m]",
                                      ytitle="$y$ [$\mu$m]",
                                      show=0,
                                      use_profiles_instead_histograms=True,
                                      xrange=[-100, 100],
                                      yrange=[-100, 100],
                                      # cmap=plt.cm.hsv,
                                      )

plt.savefig("plot_CSD_at_source.png")
plot_show()
print("File written to disk: plot_CSD_at_source.png")