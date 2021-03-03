import sys

from datastorage import DataStorage as ds
from sr import undulator
from sr import abcd
import numpy as np
from matplotlib import pyplot as plt

def source_id18(period=20, length=2.5, energy=8):
    u = undulator.get_cpmu(period=period, length=length)
    pars = u.find_harmonic_and_gap(energy)[0]
    b = u.photon_beam_characteristics(**pars)
    gsmh = abcd.GSM_Numeric(
        rms_size=b.sh,
        rms_cl=b.gsm_sclh,
        wavelen=b.wavelength * 1e-10,
    )
    gsmv = abcd.GSM_Numeric(
        rms_size=b.sv,
        rms_cl=b.gsm_sclv,
        wavelen=b.wavelength * 1e-10,
    )
    return gsmh, gsmv


def id18(
    energy=8,
    optics=[[40, "x1", "coll"], [170, "x1", "focus@200"]],
    optics_h=None,
    optics_v=None,
    samplez=200,
):

    z = np.concatenate(
        (np.arange(0, 230, 1), np.arange(samplez - 2, samplez + 2, 0.02))
    )
    z = np.unique(z)


    h, v = source_id18(period=18, length=2.5, energy=energy)
    if optics_h is None:
        optics_h = optics
    if optics_v is None:
        optics_v = optics

    v = abcd.propagate(
        beam=v,
        optics=optics_v,
        use_transfocator=False,
        z=z,
    )
        
    h = abcd.propagate(
        beam=h,
        optics=optics_h,
        use_transfocator=False,
        z=z,
    )

    ret = ds(h=h, v=v,energy=energy)
    return ret


def main(energy=8, pos_last_focusing=170):
    # each 'optics' element is
    # (distance_from_source, aperture, focal_length)
    # - aperture can be opening OR "x0.5" to indicate multiple of the coherence
    #   length at slit position
    # - focal_length can be in m or a string indicating to optimize for given condition
    r = id18(
        energy=energy,
        optics=(
            (35               , "x0.5" , None),
            (65               , None   , f"size400@{pos_last_focusing}"),
            (pos_last_focusing, None   , f"focus@200"),
        ),
    )
    # plt.plot(r.h.z,r.h.fwhm_size)
    # plt.title("horizontal")
    # print(">>>> H ", r.h.info.log)
    # plt.show()
    #
    # plt.plot(r.v.z,r.v.fwhm_size)
    # plt.title("vertical")
    # print(">>>> V ", r.v.info.log)
    # plt.show()

    from srxraylib.plot.gol import plot
    plot(
        r.h.z,r.h.fwhm_size,
        r.v.z, r.v.fwhm_size,
        legend=["H","V"],
        )

    return r

#
#########################################################################################################
#


def id18_U18(
        energy=7,
        optics=[[40, "x1", "coll"], [170, "x1", "focus@200"]],
        optics_h=None,
        optics_v=None,
        samplez=200,
):
    z = np.concatenate(
        (np.arange(0, 230, 1), np.arange(samplez - 2, samplez + 2, 0.02))
    )
    z = np.unique(z)

    h, v = source_id18(energy=energy)
    if optics_h is None:
        optics_h = optics
    if optics_v is None:
        optics_v = optics

    # beam=id18h,
    # optics=[[40, "x1", "coll"], [150, "x1", "focus@200"]],
    # z=np.arange(0, 230, 0.5),
    # use_transfocator=True,
    # transfocator=transfocator,
    # fixed_f = None,
    # fname=None,
    # force=False,

    print("\n\n\n\n\n\n")
    h = abcd.propagate(
        beam=h,
        optics=optics_h,
        use_transfocator=False,
        z=z,
    )

    print("\n\n\n\n\n\n")
    v = abcd.propagate(
        beam=v,
        optics=optics_v,
        use_transfocator=False,
        z=z,
    )



    return h, v
    # ret = ds(h=h, v=v, energy=energy)
    # return ret



def plot_loop(root="tmp"):
    # fileout = "e07keV_f2_at_170m_h.dat"
    aH = np.loadtxt("%sH.dat" % root, skiprows=1)
    aV = np.loadtxt("%sV.dat" % root, skiprows=1)
    print(aH.shape)



    f = plot(aH[:, 0], aH[:, 1],
             aV[:, 0], aV[:, 1],
             title=" trajectories", xtitle="f1 [m]", ytitle="f2 [m]",
             legend = ["H", "V"],
             show=1,)

    f = plot(aH[:, 0], aH[:, 2],
             aV[:, 0], aV[:, 2],
             aH[:, 0], aH[:, 3],
             aV[:, 0], aV[:, 3],
             title=" fwhm", xtitle="f1 [m]", ytitle="size [um]",
             legend = ["H", "V", "H at waist", "V at waist"],
             show=1,)

def plot_comparisom(root="tmp"):
    # fileout = "e07keV_f2_at_170m_h.dat"
    aH = np.loadtxt("%sH.dat" % root, skiprows=1)
    aV = np.loadtxt("%sV.dat" % root, skiprows=1)
    print(aH.shape)


    # fileout = "e07keV_f2_at_170m_h.dat"
    aWH = np.loadtxt("../ID18/e07keV_f2_at_170m_h.dat", skiprows=1)
    aWV = np.loadtxt("../ID18/e07keV_f2_at_170m_v.dat", skiprows=1)


    f = plot(aH[:, 0], aH[:, 1],
             aV[:, 0], aV[:, 1],
             aWH[:, 0], aWH[:, 1],
             aWV[:, 0], aWV[:, 1],
             title=" trajectories", xtitle="f1 [m]", ytitle="f2 [m]",
             legend = ["H GSM MARCO", "V GSM MARCO", "H WOFRY", "V WOFRY"],
             show=1,
             xrange=[0,400], yrange=[19,33])

    f = plot(aH[:, 0], aH[:, 2],
             aV[:, 0], aV[:, 2],
             # aH[:, 0], aH[:, 3],
             # aV[:, 0], aV[:, 3],
             aWH[:, 0], aWH[:, 2],
             aWV[:, 0], aWV[:, 2],
             title=" fwhm", xtitle="f1 [m]", ytitle="size [um]",
             legend = ["H GSM MARCO", "V GSM MARCO", "H WOFRY", "V WOFRY"],
             show=1,
             xrange=[0, 400], yrange=[0,75])

def run_loop(root="tmp"):

    outfileH = "%sH.dat" % root
    outfileV = "%sV.dat" % root

    if True:
        fH = open(outfileH, 'w')
        fV = open(outfileV, 'w')
        for f1 in np.linspace(10,500,200):
            h, v = id18_U18(
                energy=energy,
                optics_h=(
                    (35               , "x0.5" , None),
                    (66               , None   , f1),
                    (pos_last_focusing, None   , f"focus@200"),
                ),
                optics_v=(
                    (35, "x0.5", None),
                    (66, None, f1),
                    (pos_last_focusing, None, f"focus@200"),
                ),
            )

            # H
            i200 = np.argwhere(h.z == 200)
            iWaistH = np.argmin(h.fwhm_size)
            print("f1: %g, f2: %g" %
                  (h["info"]["focal_lengths"][2],
                  h["info"]["focal_lengths"][3]))
            print("\nSize at 200m:  H: %g" % (h.fwhm_size[i200]) )
            print("Waist H: %g (at %gm) " % (
                h.fwhm_size[iWaistH], h.z[iWaistH] ))

            fH.write("%g %g %g %g %g\n" % (
                h["info"]["focal_lengths"][2],
                h["info"]["focal_lengths"][3],
                h.fwhm_size[i200],
                h.fwhm_size[iWaistH],
                h.z[iWaistH],
            ))

            # V
            i200 = np.argwhere(v.z == 200)
            iWaistV = np.argmin(v.fwhm_size)
            print("f1: %g, f2: %g" %
                  (v["info"]["focal_lengths"][2],
                  v["info"]["focal_lengths"][3]))
            print("\nSize at 200m:  V: %g" % (v.fwhm_size[i200]) )
            print("Waist V: %g (at %gm) " % (
                v.fwhm_size[iWaistV], v.z[iWaistV] ))

            fV.write("%g %g %g %g %g\n" % (
                v["info"]["focal_lengths"][2],
                v["info"]["focal_lengths"][3],
                v.fwhm_size[i200],
                v.fwhm_size[iWaistV],
                v.z[iWaistV],
            ))

        fH.close()
        fV.close()
        print("Files written to disk: %s %s" % (outfileH, outfileV))

def run_loop_with_wofry_f1f2(root="gsm_marco_with_wofry_f1f2"):

    outfileH = "%sH.dat" % root
    outfileV = "%sV.dat" % root

    aWH = np.loadtxt("../ID18/e07keV_f2_at_170m_h.dat", skiprows=1)
    aWV = np.loadtxt("../ID18/e07keV_f2_at_170m_v.dat", skiprows=1)
    f1H = aWH[:,0]
    f2H = aWH[:,1]
    f1V = aWV[:,0]
    f2V = aWV[:,1]


    fH = open(outfileH, 'w')
    fV = open(outfileV, 'w')
    for i,tmp in enumerate(f1H):
        h, v = id18_U18(
            energy=energy,
            optics_h=(
                (35               , "x0.5" , None),
                (66               , None   , f1H[i]),
                (pos_last_focusing, None   , f2H[i]),
            ),
            optics_v=(
                (35, "x0.5", None),
                (66, None, f1V[i]),
                (pos_last_focusing, None, f2V[i]),
            ),
        )

        # H
        i200 = np.argwhere(h.z == 200)
        iWaistH = np.argmin(h.fwhm_size)
        print("f1: %g, f2: %g" %
              (h["info"]["focal_lengths"][2],
              h["info"]["focal_lengths"][3]))
        print("\nSize at 200m:  H: %g" % (h.fwhm_size[i200]) )
        print("Waist H: %g (at %gm) " % (
            h.fwhm_size[iWaistH], h.z[iWaistH] ))

        fH.write("%g %g %g %g %g\n" % (
            h["info"]["focal_lengths"][2],
            h["info"]["focal_lengths"][3],
            h.fwhm_size[i200],
            h.fwhm_size[iWaistH],
            h.z[iWaistH],
        ))

        # V
        i200 = np.argwhere(v.z == 200)
        iWaistV = np.argmin(v.fwhm_size)
        print("f1: %g, f2: %g" %
              (v["info"]["focal_lengths"][2],
              v["info"]["focal_lengths"][3]))
        print("\nSize at 200m:  V: %g" % (v.fwhm_size[i200]) )
        print("Waist V: %g (at %gm) " % (
            v.fwhm_size[iWaistV], v.z[iWaistV] ))

        fV.write("%g %g %g %g %g\n" % (
            v["info"]["focal_lengths"][2],
            v["info"]["focal_lengths"][3],
            v.fwhm_size[i200],
            v.fwhm_size[iWaistV],
            v.z[iWaistV],
        ))

    fH.close()
    fV.close()
    print("Files written to disk: %s %s" % (outfileH, outfileV))

def run_single_case():


    # each 'optics' element is
    # (distance_from_source, aperture, focal_length)
    # - aperture can be opening OR "x0.5" to indicate multiple of the coherence
    #   length at slit position
    # - focal_length can be in m or a string indicating to optimize for given condition

    # h, v = id18_U18(
    #     energy=energy,
    #     optics=(
    #         (35               , "x0.5" , None),
    #         (66               , None   , f"size400@{pos_last_focusing}"),
    #         (pos_last_focusing, None   , f"focus@200"),
    #     ),
    # )

    h, v = id18_U18(
        energy=energy,
        optics_h=(
            (35               , 36e-6 , None),
            (66               , None   , f"size400@{pos_last_focusing}"),
            (pos_last_focusing, None   , f"focus@200"),
        ),
        optics_v=(
            (35, 201e-6, None),
            (66, None, f"size400@{pos_last_focusing}"),
            (pos_last_focusing, None, f"focus@200"),
        ),
    )

    # plot(
    #     h.z, h.fwhm_size,
    #     title="H, Energy: %g keV, pos last focusing: %g m" % (energy, pos_last_focusing),
    #     )



    r = ds(h=h, v=v, energy=energy)
    #

    plot(
        r.h.z, r.h.fwhm_size,
        r.v.z, r.v.fwhm_size,
        legend=["H","V"],
        title="Energy: %g keV, pos last focusing: %g m" % (energy, pos_last_focusing),
        )

    # for i in range(len(r.h.z)):
    #     print("positions: %g, H size: %g, V size: %g" % (r.h.z[i], r.h.fwhm_size[i], r.v.fwhm_size[i]))

    # print(len(r.h.z), len(r.v.z))

    i200 = np.argwhere(r.h.z == 200)
    iWaistH = np.argmin(r.h.fwhm_size)
    iWaistV = np.argmin(r.v.fwhm_size)
    print("\nSize at 200m:  H: %g, V: %g" % (r.h.fwhm_size[i200], r.v.fwhm_size[i200]) )
    print("Waist H: %g (at %gm): V: %g (at %gm): " % (
        r.h.fwhm_size[iWaistH], r.h.z[iWaistH],
        r.v.fwhm_size[iWaistH], r.v.z[iWaistH], ))


if __name__ == "__main__":
    # main()
    from srxraylib.plot.gol import plot

    energy = 7
    pos_last_focusing = 170


    # run_single_case()
    # run_loop()
    # plot_loop()
    plot_comparisom(root="tmp")

    # run_loop_with_wofry_f1f2(root="gsm_marco_with_wofry_f1f2")
    # plot_loop(root="gsm_marco_with_wofry_f1f2")

    plot_comparisom(root="gsm_marco_with_wofry_f1f2")

