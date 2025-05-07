# -*- coding: utf-8 -*-

import numpy as np
neumann_warning = True

def fidelity(rho, phi):
    """
    Calculate the fidelity between two quantum states.
    F = <phi|rho|phi>^2 

    Parameters
    ----------
    rho : np.ndarray
        The density matrix of the first quantum state.
    phi: np.ndarray
        The reference state.

    Returns
    -------
    float
        The fidelity between the two quantum states.
    """
    F = np.dot(phi.conj().T, np.dot(rho, phi))
    return F


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

def fidelity_evolution(phi, fnpop):
        """
        Calculate the fidelity between the current and a reference state
        as a function of time.
        Parameters
        ----------
        phi : np.ndarray
            The reference state.
        fnpop : str
            The path to the population data file.
        Returns
        -------
        time : np.ndarray
            The time values from the population data file.
        F : np.ndarray
            The fidelity values between the current and reference states.
        """

        assert isinstance(phi, np.ndarray), \
            'phi must be a numpy array'
        
        time, pop = load_populations(fnpop)
        
        assert pop.shape[1:] == (phi.shape[0], phi.shape[0]), \
            'The shape of the density matrix does not match the shape of phi'
        
        # Calculate the fidelity for each time step
        F =  np.empty_like(time)
        for i, rho in enumerate(pop):
            F[i] = fidelity(rho, phi)

        return time, F


def neumann_entropy(rho):
    """
    Calculate the von Neumann entropy of a quantum state.
    S = -Tr(rho_a*log(rho_a))
    rho_a = Tr_b(rho)

    Parameters
    ----------
    rho : np.ndarray
        The density matrix of the quantum state.

    Returns
    -------
    float
        The von Neumann entropy of the quantum state.
    """
    global neumann_warning
    if neumann_warning:
        print('Warning: I don\'t know what the units of the entropy are')
        print('Warning: I truncate the density matrix to 4x4!')
        neumann_warning = False
    rho = rho[:4, :4]  
    
    #Trace
    rho_a = rho[:2, :2] + rho[2:, 2:]
    
    # Calculate the eigenvalues of the reduced density matrix
    eigvals = np.linalg.eigvalsh(rho_a)
    
    # Calculate the von Neumann entropy
    S = -np.sum(eigvals * np.log(eigvals + 1e-10))
    return S


def entropy_evolution(fnpop):
     """
    Calculate the von Neumann entropy as a function of time.
     Parameters
     ----------
     fnpop : str
         The path to the population data file.
     Returns
     -------
     time : np.ndarray
         The time values from the population data file.
     S : np.ndarray
         The von Neumann entropy values as a function of time.
     """
     time, pop = load_populations(fnpop)
     
     # Calculate the von Neumann entropy for each time step
     S = np.empty_like(time)
     for i, rho in enumerate(pop):
         S[i] = neumann_entropy(rho)
         
     return time, S