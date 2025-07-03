#!/usr/bin/env python3
# 
# Scan over combinations frequncy for different cominations of A and B 
# ./run.py H_QD-default.input Floquet-default.input

import sys
import os

import numpy as np
import matplotlib.pyplot as plt

sys.path.append('/home/matyas/Programs/PyTimeESR')
import PyTimeESR as PyTimeESR

def main(fn_ham: str, fn_dyn: str, home_path: str, code_path: str,
        Adrive: float = 0., Bdrive: float = 0., freq: float = 16.5):

    # load the input files
    ham_dict = PyTimeESR.Hamiltonian.load_input(fn_ham)
    dyn_dict = PyTimeESR.Floquet.load_input(fn_dyn)
    
    # set up parameters
    dyn_dict['n_max'] = 150
    dyn_dict['p_max'] = 210
    
    dyn_dict['A'] = [complex(0.), complex(Adrive)]
    dyn_dict['bessel_amplitude'] = [0., Bdrive]
    
    dyn_dict['frequency'] = freq
    
    # run the simulation
    Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, home_path, 
            code_path, code_version='floquet')
    Sim.run()
    
    # write output into a common file
    outfile = 'B-curr.dat'
    if not os.path.exists(outfile):
        with open(outfile, 'w') as f:
            f.write("# Adrive, Bdrive, frequency, DC\n")

    with open(outfile, 'a') as f:
        f.write(f"{Adrive}, {Bdrive}, {freq}, {Sim.results_dict['DC']}\n")

    # save output files
    oudir = os.path.join(home_path, f'freq={freq:.2f}')
    if not os.path.exists(oudir):
        os.makedirs(oudir)
    
    # what output to save
    save_output = ['current', 'populations', 'coherence', 'spin']

    for key in save_output:
        if key in Sim.output_dict:
            os.rename(Sim.output_dict[key], 
                os.path.join(oudir, f'{key}.dat'))
            

if __name__ == "__main__":
    fn_ham = sys.argv[1]
    fn_dyn = sys.argv[2]
    
    # turn into absolute path
    fn_ham = os.path.abspath(fn_ham)
    fn_dyn = os.path.abspath(fn_dyn)

    home_path = os.getcwd()
    code_path = '/home/matyas/Programs/Floquet_ESR-dev/src'

    # select amplitudes
    Bdrives = np.arange(3.,4.,.5)
    Adrives = np.arange(1.6, 2.,.1)

    # select frequencies
    freqencies = np.arange(16.5, 17.4, .05)
    freq_zoom = np.arange(16.9, 17.1, .01)
    freqencies = np.concatenate([freqencies, freq_zoom])
    freqencies = np.sort(freqencies)

    for Adrive in Adrives:
        dirA = f'A={Adrive:.3f}'
        if not os.path.exists(dirA):
            os.makedirs(dirA)

        for Bdrive in Bdrives:
            dirname = os.path.join(
                dirA, f'A={Adrive:.3f}B={Bdrive:.3f}')
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            os.chdir(dirname)
            for f in freqencies:
                main(fn_ham, fn_dyn, os.getcwd(), code_path, Adrive, Bdrive, freq=f)
            os.chdir(home_path)
