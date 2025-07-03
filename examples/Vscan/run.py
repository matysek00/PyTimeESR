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
        Vdc: float): 

    # load the input files
    ham_dict = PyTimeESR.Hamiltonian.load_input(fn_ham)
    dyn_dict = PyTimeESR.Floquet.load_input(fn_dyn)
    
    # set up parameters
    dyn_dict['bias'][1] = Vdc
    
    # run the simulation
    Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, home_path, 
            code_path, code_version='floquet')
    Sim.run()
    Sim.load_output()

    # write output into a common file
    outfile = 'Vdc-curr.dat'
    if not os.path.exists(outfile):
        with open(outfile, 'w') as f:
            f.write("# Vdc, DC\n")

    with open(outfile, 'a') as f:
        f.write(f"{Vdc}, {Sim.results_dict['DC']}\n")

            

if __name__ == "__main__":
    fn_ham = sys.argv[1]
    fn_dyn = sys.argv[2]
    
    # turn into absolute path
    fn_ham = os.path.abspath(fn_ham)
    fn_dyn = os.path.abspath(fn_dyn)

    home_path = os.getcwd()
    code_path = '/home/matyas/Programs/Floquet_ESR-dev/src'

    # select amplitudes
    Vdcs = np.arange(-15., -5., .1)

    for Vdc in Vdcs:
        main(fn_ham, fn_dyn, home_path, code_path, Vdc)