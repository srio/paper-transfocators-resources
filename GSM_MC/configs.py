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
# import transfocators


from f1_f2_scans_v3 import id18_scan,source_id18

def do_config(target=20,pos_last_focusing=170,energy=7,use_crl_lens=False):
    r=id18_scan(pos_last_focusing=pos_last_focusing, energy=energy,use_crl_lens=use_crl_lens,plot=False) 
    r = ds(r) # convert to DataStorage
    for direction in r.h,r.v: 
        idx_slope = np.gradient(direction.size_at_sample) < 0 # to avoid using short focal length with f1 focus before f2
        idx = np.argmin(np.abs(direction.size_at_sample[idx_slope]-target))
        if use_crl_lens:
            s = f"{direction.opening*1e6:.1f},{direction.f1l[idx_slope][idx]:.2f},{direction.f1_lens_radius[idx_slope][idx]*1e6:.1f},{direction.f2l[idx_slope][idx]:.2f},{direction.f2_lens_radius[idx_slope][idx]*1e6:.1f},{direction.size_at_sample[idx_slope][idx]:.1f}"
        else:
            s = f"{direction.opening*1e6:.1f},{direction.f1[idx_slope][idx]:.2f},None,{direction.f2[idx_slope][idx]:.2f},None,{direction.size_at_sample[idx_slope][idx]:.1f}"
        print(s)
        direction.out_string=s
    return r


CONFIGS = {
        1  : dict(target=20,pos_last_focusing=170,energy=7,use_crl_lens=False),
        2  : dict(target=20,pos_last_focusing=170,energy=7,use_crl_lens=True),
        3  : dict(target=30,pos_last_focusing=170,energy=7,use_crl_lens=False),
        4  : dict(target=30,pos_last_focusing=170,energy=7,use_crl_lens=True),
        5  : dict(target=7,pos_last_focusing=170,energy=7,use_crl_lens=False),
        6  : dict(target=7,pos_last_focusing=170,energy=7,use_crl_lens=True),
        7  : dict(target=12,pos_last_focusing=192,energy=7,use_crl_lens=False),
        8  : dict(target=12,pos_last_focusing=192,energy=7,use_crl_lens=True),
        9  : dict(target=3,pos_last_focusing=192,energy=7,use_crl_lens=False),
        10 : dict(target=3,pos_last_focusing=192,energy=7,use_crl_lens=True),
        }

def do_configs():
    for n,config in CONFIGS.items():
        print(config)
        do_config(**config)
        print("\n\n")



def config1(beam_source,opening=10e-6,f1=30,f2=30,pos_pinhole=36,pos_collimating=65,pos_last_focusing=170,pos_sample=200):
    beam = beam_source.propagate(pos_pinhole)
    beam = beam.gauss_aperture(opening/2.35)
    beam = beam.propagate(pos_collimating-pos_pinhole)
    beam = beam.lens(f1)
    beam = beam.propagate(pos_last_focusing-pos_collimating)
    beam = beam.lens(f2)
    beam = beam.propagate(pos_sample-pos_last_focusing)
    return beam

def do(energy=7):
    pars = dict(
            h = dict(
                E7 = dict(opening=40.0e-6, f1=41.10, f2=26.21)
            ),
            v = dict(
                E7 = dict(opening=206.8e-6, f1=49.14, f2=41.39)
            )
        )
    x = np.linspace(-4,4,401)
    for direction,beam in zip(pars.keys(),source_id18(energy=energy)):
        energy_key = f"E{energy}"
        out = config1(beam,**pars[direction][energy_key])
        plt.plot(x,out.propagate(x).rms_size*2.35*1e6,label=direction)
    plt.legend()



if __name__ == "__main__":
    #    do_all()
    do_configs()
