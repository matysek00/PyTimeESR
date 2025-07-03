#!/usr/bin/env python

import os
import sys
import numpy as np

sys.path.append('/home/matyas/Programs/PyTimeESR')

import PyTimeESR as PyTimeESR


def main():
    ham_dict = PyTimeESR.default_ham
    dyn_dict = PyTimeESR.default_dyn
    dyn_dict['runs'] = True
    
    run_path = os.getcwd()
    code_path = '/home/matyas/Programs/TimeESR/src'

    freqencies = np.arange(15., 18., 0.04)
    esr = np.zeros_like(freqencies)
    PyTimeESR.make(code_path)
    
    fnout = os.path.join(run_path, 'SpectraESR.dat')
    outfile = open(fnout, 'w')

    for i, f in enumerate(freqencies):
        dyn_dict['intervals'][0]['freq'][0]['Frequency'] = f
        
        Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, run_path, 
                                   code_path, code_version='standard')
        Sim.run()

        os.rename(
            Sim.output_dict['population'],
            os.path.join(run_path, f'POP-{f:02f}.dat'))
        os.rename(
            Sim.output_dict['spin_dynamics'],
            os.path.join(run_path, f'SP-{f:02f}.dat'))

        esr[i] = Sim.results_dict['DC']
        outfile.write(f"{f:.6f} {esr[i]:.6f}\n")
        outfile.flush()

    outfile.close()

        
if __name__ == '__main__':
    main()
