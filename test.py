#!/usr/bin/env python

import PyTimeESR as PyTimeESR
import PyTimeESR.pytimeesr

def main():
        
    ham_dict = PyTimeESR.default_ham
    dyn_dict = PyTimeESR.default_dyn
    run_path = '/home/matyas/Programs/PyTimeESR/workspace'
    code_path = '/home/matyas/Programs/TimeESR-Bessel/src'

    PyTimeESR.make(code_path)
    Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, run_path, code_path)
    Sim.run()

if __name__ == "__main__":
    main()