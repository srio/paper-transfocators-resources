#
# Import section
#
import numpy
from tools import plotCSD, plot_four_images_with_histograms, plot_show
from wofry1d_h_36m import main as main_h
from wofry1d_v_36m import main as main_v
import matplotlib.pylab as plt

#
# MAIN========================
#

from tools import plotCSD

tally_h = main_h(do_plot=0)


CSD_H, X_H, Y_H = plotCSD(
    tally_h, range_limits=[-800,800], compare_profiles=0, rotate_axes_flag=2, rotate_axes_normalization=3, normalize_to_DoC=1,
    direction='x',
    srw_file='/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/CSDx7.0keV_15k_ME_intensity.dat',
    do_plot=0)



tally_v = main_v(do_plot=0)

CSD_V, X_V, Y_V = plotCSD(
    tally_v, range_limits=[-800,800], compare_profiles=0, rotate_axes_flag=2, rotate_axes_normalization=3, normalize_to_DoC=1,
    direction='y',
    srw_file='/users/srio/OASYS1.2/paper-transfocators-resources/workspaces/SRW/CSD/CSDy7.0keV_15k_ME_intensity.dat',
    do_plot=0)

CSD_H.append(CSD_V[0])
CSD_H.append(CSD_V[1])
X_H.append(X_V[0])
X_H.append(X_V[1])
Y_H.append(Y_V[0])
Y_H.append(Y_V[1])

plot_four_images_with_histograms(CSD_H, X_H, Y_H,
                  title_list=["WOFRY1D H", "SRW-ME H", "WOFRY1D V", "SRW-ME V"],
                  xtitle=["$(x_1+x_2)/2$ [$\mu$m]", "$(x_1+x_2)/2$ [$\mu$m]", "$(y_1+y_2)/2$ [$\mu$m]", "$(y_1+y_2)/2$ [$\mu$m]",],
                  ytitle=["$(x_2-x_1)$ [$\mu$m]", "$(x_2-x_1)$ [$\mu$m]", "$(y_2-y_1)$ [$\mu$m]", "$(y_2-y_1)$ [$\mu$m]",],
                  show=0,
                  use_profiles_instead_histograms=True,
                  xrange=[-800, 800],
                  yrange=[-800, 800],
                  # cmap=plt.cm.hsv,
                  )

plt.savefig("plot_DoC_at_36m.png")
plot_show()
print("File written to disk: plot_DoC_at_36m.png")