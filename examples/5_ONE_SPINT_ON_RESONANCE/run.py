#!/usr/bin/env python

import os
import sys
import numpy as np

sys.path.append('/home/matyas/Programs/PyTimeESR')

import PyTimeESR as PyTimeESR


def main():
    ham_dict = PyTimeESR.default_ham
    dyn_dict = PyTimeESR.default_dyn
    
    ham_dict['eps_QD'] = -5.0
    ham_dict['N_plot'] = 100

    dyn_dict['Ntime'] = 150000
    dyn_dict['intervals'][0]['freq'][0]['Frequency'] = 17.02717787
    dyn_dict['gamma0'] = [.001, .0005]
    dyn_dict['gamma1'] = [.0, .0001]
    dyn_dict['gammaC'] = 0.001

    dyn_dict['biases'][0]['bias'] = [3., -3.]
    dyn_dict['Temperature'] = 0.05
    dyn_dict['density_matrix'] = True
    dyn_dict['redimension'] = True
    dyn_dict['Nd'] = 3
    
    run_path = os.getcwd()
    code_path = '/home/matyas/Programs/TimeESR-Bessel/src'

    
    PyTimeESR.make(code_path)
        
    Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, run_path, code_path)
    Sim.run()

        
if __name__ == '__main__':
    main()
