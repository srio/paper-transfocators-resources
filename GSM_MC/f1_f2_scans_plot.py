import sys

from datastorage import DataStorage as ds
from sr import undulator
from sr import abcd
import numpy as np
from matplotlib import pyplot as plt


def plot_comparison_f1f2_scan(root="f1_f2_scan_7keV_f2_at_170", f1f2_from='each'):

    if f1f2_from == "each":
        aH = np.loadtxt("%s_h.dat" % root, skiprows=1)
        aV = np.loadtxt("%s_v.dat" % root, skiprows=1)

        aWH = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_h.dat", skiprows=1)
        aWV = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_v.dat", skiprows=1)

    elif f1f2_from == "marco":
        aH = np.loadtxt("%s_h.dat" % root, skiprows=1)
        aV = np.loadtxt("%s_v.dat" % root, skiprows=1)

        aWH = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_h_f1f2_from_marco.dat", skiprows=1)
        aWV = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_v_f1f2_from_marco.dat", skiprows=1)

    elif f1f2_from == "wofry":
        aH = np.loadtxt("%s_h_f1f2_from_wofry.dat" % root, skiprows=1)
        aV = np.loadtxt("%s_v_f1f2_from_wofry.dat" % root, skiprows=1)

        aWH = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_h.dat", skiprows=1)
        aWV = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_v.dat", skiprows=1)



    f = plot(aH[:, 0], aH[:, 1],
             aV[:, 0], aV[:, 1],
             aWH[:, 0], aWH[:, 1],
             aWV[:, 0], aWV[:, 1],
             title=" trajectories from "+f1f2_from, xtitle="f1 [m]", ytitle="f2 [m]",
             legend = ["H GSM MARCO", "V GSM MARCO", "H WOFRY", "V WOFRY"],
             show=1,
             xrange=[0,100], yrange=[5,100])

    f = plot(aH[:, 0], aH[:, 2],
             aV[:, 0], aV[:, 2],
             # aH[:, 0], aH[:, 3],
             # aV[:, 0], aV[:, 3],
             aWH[:, 0], aWH[:, 2],
             aWV[:, 0], aWV[:, 2],
             title=" trajectories from "+f1f2_from, xtitle="f1 [m]", ytitle="size FWHM [um]",
             legend = ["H GSM MARCO", "V GSM MARCO", "H WOFRY", "V WOFRY"],
             show=1,
             xrange=[0, 100], yrange=[0,75])




if __name__ == "__main__":
    # main()
    from srxraylib.plot.gol import plot, set_qt
    set_qt()

    energy = 7
    pos_last_focusing = 170

    plot_comparison_f1f2_scan(root="f1_f2_scan_7keV_f2_at_170", f1f2_from='each')
    # plot_comparison_f1f2_scan(root="f1_f2_scan_7keV_f2_at_170", f1f2_from='wofry')
    # plot_comparison_f1f2_scan(root="f1_f2_scan_7keV_f2_at_170", f1f2_from='marco')