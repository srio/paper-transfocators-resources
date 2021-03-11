import sys

import os
from datastorage import DataStorage as ds
import datastorage
from sr import undulator
from sr import abcd
from sr.abcd.propagate import find_fl_to_get_size
import numpy as np
from matplotlib import pyplot as plt
import scipy

def source_id18(period=20, length=2.5, energy=8, fully_coherent=False):
    u = undulator.get_cpmu(period=period, length=length,min_gap=4.8)
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
    if fully_coherent:
        gsmh = gsmh.coherent()
        gsmv = gsmv.coherent()
    return gsmh, gsmv


def print_ds(d):
    for key, value in d.items():
        print(f"{key:.20s} : {str(value)}")



def id18_scan(
    energy=7,
    pos_pinhole=35,
    pinhole=0.5, # in unit of FWHM coherence length
    pos_collimating=66,
    pos_last_focusing=170,
    pos_sample=200,
    use_gaussian_aperture=True,
    fname=None,
    plot=True,
):
    source = source_id18(energy=energy)

    f1_list = np.linspace(2,100,400)

    out = dict(h=dict(f1=[],f2=[],size_at_f2=[],size_at_sample=[],slit=None),
               v=dict(f1=[],f2=[],size_at_f2=[],size_at_sample=[],slit=None))
    for beam,store in zip(source,out.values()):
        b1 = beam.propagate(pos_pinhole)
        opening = b1.rms_cl*2.35*pinhole
        if use_gaussian_aperture:
            gaussian_opening = opening/2.35 # following Manuel, gaussian fwhm = hard opening
            aperture = abcd.GaussianAperture(gaussian_opening)
        else:
            raise ValueError("HardAperture to be implemented following Manuel CF matching")
            aperture = abcd.HardAperture(gaussian_opening)
        b2 = b1.apply(aperture)
        b3 = b2.propagate(pos_collimating-pos_pinhole)
        store["opening"] = opening
        for f1 in f1_list:
            b4 = b3.lens(f1)
            b5 = b4.propagate(pos_last_focusing-pos_collimating)
            f2 = find_fl_to_get_size(b5,pos_sample-pos_last_focusing,0,verbose=False)
            b6 = b5.lens(f2).propagate(pos_sample-pos_last_focusing)
            store["f1"].append(f1)
            store["f2"].append(f2)
            store["size_at_sample"].append(b6.rms_size*2.35*1e6)
            store["size_at_f2"].append(b5.rms_size*2.35*1e6)

    # done calculating, now work on plot
    if plot:
        fig,ax = plt.subplots(3,1,figsize=[6,8],sharex=True)
        for d in "h","v":
            ax[0].plot(out[d]["f1"],out[d]["f2"],label=d)
            ax[1].plot(out[d]["f1"],out[d]["size_at_f2"],label=d)
            ax[2].plot(out[d]["f1"],out[d]["size_at_sample"],label=d)
        ax[2].set_xlabel("f1 [m]")
        ax[0].set_ylabel("f2 [m]")
        ax[1].set_ylabel("fwhm size at f2 [um]")
        ax[2].set_ylabel("fwhm size at sample [um]")
        for a in ax: a.grid()
        sh = out['h']["opening"]*1e6*2.35
        sv = out['v']["opening"]*1e6*2.35
        title = f"E={energy}, pos (slit,f1,f2)=({pos_pinhole:.0f},{pos_collimating:.0f},{pos_last_focusing:.0f}), slit = {sh:.0f},{sv:.0f}"
        ax[0].set_title(title)
        ax[1].set_ylim(0,500)
        if pos_last_focusing == 170:
            ax[0].set_ylim(15,40)
        else:
            ax[0].set_ylim(5,20)
        plt.tight_layout()

    # write to file (added srio)

    for direction in ['h','v']:

        filename = 'f1_f2_scan_%dkeV_f2_at_%d_%s.dat' % (energy, pos_last_focusing, direction)
        # filename = 'f1_f2_scan_%dkeV_f2_at_%d_%s_f1f2_from_wofry.dat' % (energy, pos_last_focusing, direction)
        f = open(filename, 'w')
        for i in range(len(out['h']['f1'])):
            f.write("%g %g %g %g %g\n" % (
                                       out[direction]['f1'][i],
                                       out[direction]['f2'][i],
                                       out[direction]['size_at_sample'][i],
                                       out[direction]['size_at_f2'][i],
                                       out[direction]['opening'],
            ))
        f.close()
        print("File written to disk: %s " % filename)

    return out



def id18_scan_from_wofry(
    energy=7,
    pos_pinhole=35,
    pinhole=0.5, # in unit of FWHM coherence length
    pos_collimating=66,
    pos_last_focusing=170,
    pos_sample=200,
    use_gaussian_aperture=True,
    fname=None,
    plot=True,
):
    source = source_id18(energy=energy)


    a_h = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_h.dat", skiprows=1)
    a_v = np.loadtxt("../ID18_U18/e07keV_f2_at_170m_v.dat", skiprows=1)

    f1_list_h = a_h[:, 0]
    f1_list_v = a_v[:, 0]
    f2_list_h = a_h[:, 1]
    f2_list_v = a_v[:, 1]
    F1_list = [f1_list_h, f1_list_v]
    F2_list = [f2_list_h, f2_list_v]


    # f1_list = np.linspace(2,100,400)

    out = dict(h=dict(f1=[],f2=[],size_at_f2=[],size_at_sample=[],slit=None),
               v=dict(f1=[],f2=[],size_at_f2=[],size_at_sample=[],slit=None))

    for beam,store,f1_list,f2_list in zip(source,out.values(),F1_list,F2_list):
    # for beam,store in zip(source,out.values()):
        b1 = beam.propagate(pos_pinhole)
        opening = b1.rms_cl*2.35*pinhole
        if use_gaussian_aperture:
            gaussian_opening = opening/2.35 # following Manuel, gaussian fwhm = hard opening
            aperture = abcd.GaussianAperture(gaussian_opening)
        else:
            raise ValueError("HardAperture to be implemented following Manuel CF matching")
            aperture = abcd.HardAperture(gaussian_opening)
        b2 = b1.apply(aperture)
        b3 = b2.propagate(pos_collimating-pos_pinhole)
        store["opening"] = opening
        for f1,f2 in zip(f1_list,f2_list):
            b4 = b3.lens(f1)
            b5 = b4.propagate(pos_last_focusing-pos_collimating)
            # f2 = find_fl_to_get_size(b5,pos_sample-pos_last_focusing,0,verbose=False)
            b6 = b5.lens(f2).propagate(pos_sample-pos_last_focusing)
            store["f1"].append(f1)
            store["f2"].append(f2)
            store["size_at_sample"].append(b6.rms_size*2.35*1e6)
            store["size_at_f2"].append(b5.rms_size*2.35*1e6)

    # done calculating, now work on plot
    if plot:
        fig,ax = plt.subplots(3,1,figsize=[6,8],sharex=True)
        for d in "h","v":
            ax[0].plot(out[d]["f1"],out[d]["f2"],label=d)
            ax[1].plot(out[d]["f1"],out[d]["size_at_f2"],label=d)
            ax[2].plot(out[d]["f1"],out[d]["size_at_sample"],label=d)
        ax[2].set_xlabel("f1 [m]")
        ax[0].set_ylabel("f2 [m]")
        ax[1].set_ylabel("fwhm size at f2 [um]")
        ax[2].set_ylabel("fwhm size at sample [um]")
        for a in ax: a.grid()
        sh = out['h']["opening"]*1e6*2.35
        sv = out['v']["opening"]*1e6*2.35
        title = f"E={energy}, pos (slit,f1,f2)=({pos_pinhole:.0f},{pos_collimating:.0f},{pos_last_focusing:.0f}), slit = {sh:.0f},{sv:.0f}"
        ax[0].set_title(title)
        ax[1].set_ylim(0,500)
        if pos_last_focusing == 170:
            ax[0].set_ylim(15,40)
        else:
            ax[0].set_ylim(5,20)
        plt.tight_layout()

    # write to file (added srio)

    for direction in ['h','v']:
        filename = 'f1_f2_scan_%dkeV_f2_at_%d_%s_f1f2_from_wofry.dat' % (energy, pos_last_focusing, direction)
        f = open(filename, 'w')
        for i in range(len(out['h']['f1'])):
            f.write("%g %g %g %g %g\n" % (
                                       out[direction]['f1'][i],
                                       out[direction]['f2'][i],
                                       out[direction]['size_at_sample'][i],
                                       out[direction]['size_at_f2'][i],
                                       out[direction]['opening'],
            ))
        f.close()
        print("File written to disk: %s " % filename)

    return out



def do_all(energy=7):
    # o=id18_scan(pos_last_focusing=190,energy=energy)
    # plt.savefig("/tmp/f1_f2_scan_190m.png",transparent=True,dpi=300)

    # o=id18_scan(pos_last_focusing=170,energy=energy)
    # plt.savefig("/tmp/f1_f2_scan_170m.png",transparent=True,dpi=300)
    # plt.show()

    o=id18_scan_from_wofry(pos_last_focusing=170,energy=energy)
    # plt.savefig("/tmp/f1_f2_scan_170m.png",transparent=True,dpi=300)
    plt.show()


if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()

    do_all()
    # pass
