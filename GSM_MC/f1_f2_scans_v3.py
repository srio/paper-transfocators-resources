import sys

import os
from datastorage import DataStorage as ds
import datastorage
import sr
from sr import undulator
from sr import abcd
from sr.abcd.propagate import find_fl_to_get_size
import numpy as np
from matplotlib import pyplot as plt
import scipy
import transfocators


def source_id18(period=18, length=2.5, energy=8, fully_coherent=False):
    u = undulator.get_cpmu(period=period, length=length, min_gap=4.8)
    pars = u.find_harmonic_and_gap(energy,sort_harmonics=True)[0]
    b = u.photon_beam_characteristics(**pars)
    gsmh = abcd.GSM_Numeric(
        rms_size=b.sh, rms_cl=b.gsm_sclh, wavelen=b.wavelength * 1e-10,
    )
    gsmv = abcd.GSM_Numeric(
        rms_size=b.sv, rms_cl=b.gsm_sclv, wavelen=b.wavelength * 1e-10,
    )
    if fully_coherent:
        gsmh = gsmh.coherent()
        gsmv = gsmv.coherent()
    return gsmh, gsmv


def print_ds(d):
    for key, value in d.items():
        print(f"{key:.20s} : {str(value)}")


def get_lens(f, energy, use_crl_aperture=True):
    """ return sr.crl.LensBlock with
    radius matching the focal lenght
    """
    delta, beta = sr.crl.get_delta_beta("Be", energy=energy)
    radius = f * 2 * delta
    lens = sr.crl.LensBlock(
        n=1,
        radius=radius,
        material="Be",
        thickness=0.001,
        #        thickness=0.01,
        web_thickness=3e-05,
    )
    if not use_crl_aperture:
        lens._aperture = 1e3
    return lens


def add_transfocator_plot(ax, title, energy, out, pos_collimating, pos_last_focusing):
    out = ds(out)  # symplify typing out["h"]["f1"]→ out.h.f1
    t1 = transfocators.TRANSFOCATORS_2D[pos_collimating]
    t2 = transfocators.TRANSFOCATORS_2D[pos_last_focusing]
    t1h = transfocators.TRANSFOCATORS_H[pos_collimating]
    t2h = transfocators.TRANSFOCATORS_H[pos_last_focusing]
    color_h = ax[0].lines[0].get_color()
    color_v = ax[0].lines[1].get_color()

    df1f2 = out.h.df1f2
    df2sample = out.h.df2sample

    newfig,newax=plt.subplots(1,1)

    # find range of horizontal beam size (h is less tunable)
    hmin = out.h.size_at_sample[-1]
    hmax = np.max(out.h.size_at_sample) * 0.9

    targets = np.linspace(hmin, hmax, 5)
    colors = plt.rcParams["axes.prop_cycle"][2:]
    for target, color in zip(targets, colors):
        color = color["color"]
        idx = np.argwhere(out.v.size_at_sample > target).ravel()[-1] + 1
        if idx >= len(out.v.f1):
            idx = len(out.v.f1) - 1
        f1_2d = out.v.f1[idx]
        f2_2d = out.v.f2[idx]
        ax[2].plot(f1_2d, target, "o", color=color)
        f1vt = t1.find_best_set_for_focal_length(energy, f1_2d)
        f2vt = t2.find_best_set_for_focal_length(energy, f2_2d)

        idx = np.argwhere(out.h.size_at_sample > target).ravel()[-1] + 1
        f1h = out.h.f1[idx]
        f2h = out.h.f2[idx]

        f1_1D = 1 / (1 / f1h - 1 / f1vt.focal_length)
        if f1_1D < 0: f1_1D = np.inf
        f2_1D = 1 / (1 / f2h - 1 / f2vt.focal_length)
        if f2_1D < 0: f2_1D = np.inf
        # print(f1_1D,f2_1D)
        f1_1Dt = t1h.find_best_set_for_focal_length(energy, f1_1D)
        f2_1Dt = t2h.find_best_set_for_focal_length(energy, f2_1D)

        print("f1",str(f1_1Dt.best_lens_set))
        print("f2",str(f2_1Dt.best_lens_set))
        print(f"f1_1D target,value = {f1_1D:.2f},{f1_1Dt.focal_length:.2f}")
        print(f"f2_1D target,value = {f2_1D:.2f},{f2_1Dt.focal_length:.2f}")
        # add again 2D lens for plotting porpouse
        f1ht = 1 / (1 / f1_1Dt.focal_length + 1 / f1vt.focal_length)
        f2ht = 1 / (1 / f2_1Dt.focal_length + 1 / f2vt.focal_length)
        print(f"f1_h target,value = {f1h:.2f},{f1ht:.2f}")
        print(f"f2_h target,value = {f2h:.2f},{f2ht:.2f}")
        ax[0].plot(f1ht, f2ht, ".", color=color)
        ax[0].plot(f1vt.focal_length, f2vt.focal_length, ".", color=color)
        print(target, f1_2d, f2_2d)
        sv = out.v.b3.apply(f1vt.best_lens_set).propagate(df1f2).apply(f2vt.best_lens_set).propagate(df2sample).rms_size*2.35*1e6
        sh = out.h.b3.apply(f1vt.best_lens_set).apply(f1_1Dt.best_lens_set).propagate(df1f2).apply(f2vt.best_lens_set).apply(f2_1Dt.best_lens_set).propagate(df2sample).rms_size*2.35*1e6
        newax.plot(target,sv,"o",color=color_v)
        newax.plot(target,sh,"o",color=color_h)
    x = np.linspace(0,hmax*1.2,100)
    newax.grid()
    newax.set_title(title)
    newax.set_xlabel("requested beamsize [μm]")
    newax.set_ylabel("obtained beamsize [μm]")
    newax.legend(["v","h"])
    newax.plot(x,x,color="0.5")
    newfig.tight_layout()
    return newfig



def add_transfocator_plot_old(ax, energy, out, pos_collimating, pos_last_focusing):
    out = ds(out)  # symplify typing out["h"]["f1"]→ out.h.f1
    t1 = transfocators.TRANSFOCATORS_2D[pos_collimating]
    t2 = transfocators.TRANSFOCATORS_2D[pos_last_focusing]
    t1h = t1 + transfocators.TRANSFOCATORS_H[pos_collimating]
    t2h = t2 + transfocators.TRANSFOCATORS_H[pos_last_focusing]
    color_h = ax.lines[0].get_color()
    color_v = ax.lines[1].get_color()
    for f1v, f2v, f1h, f2h in zip(out.v.f1, out.v.f2, out.h.f1, out.h.f2):
        f1vt = t1.find_best_set_for_focal_length(energy, f1v).focal_length
        f2vt = t2.find_best_set_for_focal_length(energy, f2v).focal_length
        ax.plot(f1vt, f2vt, ".", color=color_v)
        f1ht = t1h.find_best_set_for_focal_length(energy, f1h).focal_length
        f2ht = t2h.find_best_set_for_focal_length(energy, f2h).focal_length
        ax.plot(f1ht, f2ht, ".", color=color_h)


def id18_scan(
    energy=7,
    pos_pinhole=36,
    pinhole=0.5,  # in unit of FWHM coherence length
    pos_collimating=65,
    pos_last_focusing=170,
    pos_sample=200,
    use_gaussian_aperture=True,
    use_crl_lens=True,
    use_crl_aperture=True,
    fname=None,
    plot=True,
    plot_transfocator=True,
    specfileroot="",
):
    source = source_id18(energy=energy)

    f1_list = np.linspace(2, 100, 1001)
    n = len(f1_list)

    out = dict(
        h=dict(
            f1=np.zeros(n),
            f2=np.zeros(n),
            size_at_f2=np.zeros(n),
            size_at_sample=np.zeros(n),
            slit=None,
            b3=None,
            df1f2=pos_last_focusing - pos_collimating,
            df2sample=pos_sample - pos_last_focusing,
        ),
        v=dict(
            f1=np.zeros(n),
            f2=np.zeros(n),
            size_at_f2=np.zeros(n),
            size_at_sample=np.zeros(n),
            slit=None,
            b3=False,
            df1f2=pos_last_focusing - pos_collimating,
            df2sample=pos_sample - pos_last_focusing,
        ),
    )
    for beam, store in zip(source, out.values()):
        b1 = beam.propagate(pos_pinhole)
        opening = b1.rms_cl * 2.35 * pinhole
        if use_gaussian_aperture:
            gaussian_opening = (
                opening / 2.35
            )  # following Manuel, gaussian fwhm = hard opening
            aperture = abcd.GaussianAperture(gaussian_opening)
        else:
            raise ValueError(
                "HardAperture to be implemented following Manuel CF matching"
            )
            aperture = abcd.HardAperture(gaussian_opening)
        b2 = b1.apply(aperture)
        b3 = b2.propagate(pos_collimating - pos_pinhole)
        store["b3"]=b3
        store["opening"] = opening
        for if1,f1 in enumerate(f1_list):
            if use_crl_lens:
                f1l = get_lens(f1, energy, use_crl_aperture=use_crl_aperture)
                b4 = b3.apply(f1l)
            else:
                b4 = b3.lens(f1)
            b5 = b4.propagate(pos_last_focusing - pos_collimating)
            f2 = find_fl_to_get_size(
                b5, pos_sample - pos_last_focusing, 0, verbose=False
            )
            if use_crl_lens:
                f2l = get_lens(f2, energy, use_crl_aperture=use_crl_aperture)
                b6 = b5.apply(f2l).propagate(pos_sample - pos_last_focusing)
            else:
                b6 = b5.lens(f2).propagate(pos_sample - pos_last_focusing)
            store["f1"][if1]=f1
            store["f2"][if1]=f2
            store["size_at_sample"][if1] = b6.rms_size * 2.35 * 1e6
            store["size_at_f2"][if1] = b5.rms_size * 2.35 * 1e6


    #
    # srio  dump
    #


    if specfileroot != "":
        for d in "h", "v":
            f1 = out[d]["f1"]
            f2 = out[d]["f2"]
            size_at_f2 = out[d]["size_at_f2"]
            size_at_sample = out[d]["size_at_sample"]


            f = open("%s_%s.dat" % (specfileroot, d),'w')
            f.write("#S data\n#N 4\n#L f1  f2  size_at_f2  size_at_sample\n")
            for i in range(f1.size):
                f.write("%g  %g  %g  %g \n" % (f1[i], f2[i], size_at_f2[i], size_at_sample[i]))
            f.close()
            print("File written to disK: ", "%s_%s.dat" % (specfileroot, d))


    # done calculating, now work on plot
    if plot:
        fig, ax = plt.subplots(3, 1, figsize=[6, 8], sharex=True)
        for d in "h", "v":
            ax[0].plot(out[d]["f1"], out[d]["f2"], label=d)
            ax[1].plot(out[d]["f1"], out[d]["size_at_f2"], label=d)
            ax[2].plot(out[d]["f1"], out[d]["size_at_sample"], label=d)
        ax[2].set_xlabel("f1 [m]")
        ax[0].set_ylabel("f2 [m]")
        ax[1].set_ylabel("fwhm size at f2 [um]")
        ax[2].set_ylabel("fwhm size at sample [um]")
        for a in ax:
            a.grid()
        sh = out["h"]["opening"] * 1e6 # * 2.35
        sv = out["v"]["opening"] * 1e6 #* 2.35
        title = f"E={energy}, pos (slit,f1,f2)=({pos_pinhole:.0f},{pos_collimating:.0f},{pos_last_focusing:.0f}), slit = {sh:.0f},{sv:.0f}"
        ax[0].set_title(title)
        ax[1].set_ylim(0, 500)
        if pos_last_focusing == 170:
            ax[0].set_ylim(15, 40)
        else:
            ax[0].set_ylim(6.5, 9)
        ax[2].set_ylim(0, None)
        ax[2].legend()
        if plot_transfocator and use_crl_lens:
            fig2=add_transfocator_plot(ax, title,energy, out, pos_collimating, pos_last_focusing)
            fig2.savefig(f"results/results_{energy:.1f}keV_f2_at_{pos_last_focusing}m_beamsize_scan.png",transparent=True,dpi=300)
        fig.tight_layout()
        fig.savefig(f"results/results_{energy:.1f}keV_f2_at_{pos_last_focusing}m_f1_f2_scan.png",transparent=True,dpi=300)
    return out


def compare_apertures(energy=7, pos_last_focusing=190):
    o1 = id18_scan(
        pos_last_focusing=pos_last_focusing,
        energy=energy,
        plot=False,
        use_crl_lens=False,
    )
    o2 = id18_scan(
        pos_last_focusing=pos_last_focusing,
        energy=energy,
        plot=False,
        use_crl_lens=True,
        use_crl_aperture=False,
    )
    o3 = id18_scan(
        pos_last_focusing=pos_last_focusing,
        energy=energy,
        plot=False,
        use_crl_lens=True,
        use_crl_aperture=True,
    )
    fig, ax = plt.subplots(2, 1, figsize=[6, 8], sharex=True)
    labels = (
        "no abs, no finite aperture",
        "   abs, no finite aperture",
        "   abs,    finite aperture",
    )

    for o, l in zip((o1, o2, o3), labels):
        ax[0].plot(o["h"]["f1"], o["h"]["size_at_sample"], label=l)
        ax[1].plot(o["v"]["f1"], o["v"]["size_at_sample"], label=l)
    ax[0].legend()
    ax[1].set_xlabel("f1 [m]")
    ax[0].set_ylabel("fwhm size at sample H [um]")
    ax[1].set_ylabel("fwhm size at sample V [um]")
    for a in ax:
        a.grid()
    title = f"E={energy}, pos f2={pos_last_focusing:.0f}"
    ax[0].set_title(title)
    plt.tight_layout()
    fname = f"results/{energy:.1f}keV_f2_at_{pos_last_focusing}m_lens_aperture.png"
    plt.savefig(fname, transparent=True, dpi=300)


def do_f1_f2_scans(energy=7):
    o = id18_scan(pos_last_focusing=190, energy=energy)
    plt.savefig("/tmp/f1_f2_scan_190m.png", transparent=True, dpi=300)
    o = id18_scan(pos_last_focusing=170, energy=energy)
    plt.savefig("/tmp/f1_f2_scan_170m.png", transparent=True, dpi=300)


def do_all(close_plot=True):
    for energy in (7,10,15,20,25,30,35):
        for dist in (170,192):
            o = id18_scan(pos_last_focusing=dist, energy=energy)
            compare_apertures(energy=energy,pos_last_focusing=dist)
        if close_plot: plt.close("all")


if __name__ == "__main__":
    #    do_all()
    # pass

    o = id18_scan(pos_last_focusing=170, energy=7, specfileroot="f1_f2_scans_170")
    o = id18_scan(pos_last_focusing=192, energy=7, specfileroot="f1_f2_scans_192")
    # compare_apertures(energy=7, pos_last_focusing=170)
    # plt.close("all")
