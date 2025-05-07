# -*- coding: utf-8 -*-

import os
import numpy as np

def make(path): 
    """Compile the TimeESR code.
    """
    pwd = os.getcwd()
    os.chdir(path)
    os.system('make clean')
    os.system('make')
    os.chdir(pwd)


def load_populations(fnpop):
    """
    Load the populations from a file.
    
    Parameters
    ----------
    fnpop : str
        The path to the population data file.
    
    Returns
    -------
    time : np.ndarray
        The time values from the population data file.
    pop : np.ndarray
        The populations from the population data file.
    """
    
    pop = np.loadtxt(fnpop)
    n_states = int(np.sqrt((pop.shape[1] - 1)/2))
        
    assert 2*(n_states**2) + 1 == pop.shape[1], \
            'The population data file does not have the correct number of columns'
        
    time = pop[:, 0]

        # Extract the real and imaginary parts of the density matrix
    real_pop = pop[:, 1:n_states**2 + 1]
    imag_pop = pop[:, n_states**2 + 1:]
    pop = real_pop + 1j*imag_pop
        
        # and reshape it to the correct size
    pop = pop.reshape((-1, n_states, n_states))
    
    return time, pop
