# -*- coding: utf-8 -*-

import numpy as np
from .utils import load_populations

neumann_warning = True
concurence_warning = True
Pauli_Matrices = np.array([[[0, 1], [1, 0]],
                           [[0, -1j], [1j, 0]],
                            [[1, 0], [0, -1]],
                           ], dtype=complex)

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

def Concurrence(rho):
    """
    Calculate the concurrence of a quantum state.
    C = max(0, sqrt(lambda_1) - sqrt(lambda_2) - sqrt(lambda_3) - sqrt(lambda_4))

    Parameters
    ----------
    rho : np.ndarray
        The density matrix of the quantum state.

    Returns
    -------
    float
        The concurrence of the quantum state.
    """
    
    global concurence_warning
    if concurence_warning:
        print('Warning: I truncate the density matrix to 4x4!')
        concurence_warning = False
    rho = rho[:4, :4]
    
    #yy_mat =np.fliplr(np.diag([-1, 1, 1, -1]))
    #sigma = rho.dot(yy_mat).dot(rho.conj()).dot(yy_mat)
    #w = np.sort(np.real(np.linalg.eigvals(sigma)))
    #w = np.sqrt(np.maximum(w, 0.))
    #return max(0.0, w[-1] - np.sum(w[0:-1]))

    # Flip the density matrix spin
    transform = np.fliplr(np.diag([-1, 1, 1, -1]))#np.outer(Pauli_Matrices[1], Pauli_Matrices[1])
    rho_flipped = transform.dot(rho.conj()).dot(transform)
        
    A = rho.dot(rho_flipped)
    

    # Calculate the eigenvalues 
    eigvals = np.linalg.eigvals(A)
    eigvals = np.sort(np.sqrt(np.maximum(np.real(eigvals),0.)))
    

    # Calculate the concurrence
    C = max(0., eigvals[3] - eigvals[:3].sum())
    
    return C


def Entanglement(C):
    """
    Calculate the entanglement of a quantum state.
    E = -C*log2(C) - (1-C)*log2(1-C)
    Parameters
    ----------
    C : float
        The concurrence of the quantum state.
    Returns
    -------
    float
        The entanglement of the quantum state.
    """

    E = -C*np.log2(C) - (1-C)*np.log2(1-C)
    return E


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

def time_evolution(fn: callable, fnpop: str, *args, **kwargs):
    """
    Calculate time evolution for a given function.
    Parameters
    ----------
    fn : function
        The function to be evaluated.
    fnpop : str
        The path to the population data file.
    args : tuple
        The arguments to be passed to the function.
    kwargs : dict
        The keyword arguments to be passed to the function.
    Returns
    -------
    time : np.ndarray
        The time values from the population data file.
    result : np.ndarray
        The result of the function evaluated at each time step.
    """
    time, pop = load_populations(fnpop)
    
    # Calculate the result for each time step
    result = np.empty_like(time)
    for i, rho in enumerate(pop):
        result[i] = fn(rho, *args, **kwargs)
        
    return time, result