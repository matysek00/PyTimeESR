#!/usr/bin/env python

import os
import sys
import numpy as np

sys.path.append('/home/matyas/Programs/PyTimeESR')

import PyTimeESR


def main():
    ham_dict = PyTimeESR.default_ham
    dyn_dict = PyTimeESR.default_dyn
    home_path = os.getcwd()
    code_path = '/home/matyas/Programs/TimeESR-Bessel/src'

    dyn_dict['Ntime'] = int(10e5)
    
    dyn_dict['gamma0'] = [0.005, 0.0001]
    dyn_dict['gamma1'] = [0.0000, 0.00015]
    dyn_dict['gammaC'] = 0.005
    dyn_dict['runs'] = True
    
    freqstart = 16.8
    freqstop = 17.2
    freqdiff = [.02, .04, .08, .16]
    
    PyTimeESR.make(code_path)
    
    for fd in freqdiff: 
        freqencies = np.arange(freqstart, freqstop, fd)
        esr = np.zeros_like(freqencies)
        run_path = os.path.join(home_path, f'fd-{fd:02f}')
        
        os.mkdir(run_path)
        os.chdir(run_path)

        fnout = os.path.join(run_path, 'SpectraESR.dat')
        outfile = open(fnout, 'w')
        
        for i, f in enumerate(freqencies):
            dyn_dict['intervals'][0]['freq'][0]['Frequency'] = f
            
            Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, run_path, code_path)
            Sim.run()

            os.rename(
                Sim.output_dict['population'],
                os.path.join(run_path, f'POP-{f:02f}.dat')
            )
            os.rename(
                Sim.output_dict['spin_dynamics'],
                os.path.join(run_path, f'SP-{f:02f}.dat')
            )

            esr[i] = np.loadtxt(
                Sim.output_dict['esr']
            )
            outfile.write(f"{f:.6f} {esr[i]:.6f}\n")

        os.chdir(home_path)
        
if __name__ == '__main__':
    main()
