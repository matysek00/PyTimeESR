import numpy as np
from scipy.optimize import curve_fit



def create_fit_esr(x, y, beta0=None, lb=None, ub=None):
    """Probably redundant funcition
    """
    # Perform curve fitting
    popt, pcov = curve_fit(ESR_fit_fun, x, y, p0=beta0, bounds=(lb, ub))   
    
    # Return outputs
    return popt, np.sqrt(np.diag(pcov))

def ESR_fit_fun(x, res, gamma, Isym, Iasym, p0, p1, p2, pinv):
    """
    ESR fitting function.

    Parameters:
    x : array_like
        The independent variable (frequency or field).
    res : float
        Resonance position.
    gamma : float
        Lorentzian width (damping).
    Isym : float
        Symmetric peak current amplitude.
    Iasym : float
        Asymmetric peak current amplitude (Fano).
    p0, p1, p2 : float
        Coefficients for the polynomial background current.
    pinv : float
        Inverse background current coefficient.
    Returns:
    y : array_like
        The fitted values for the given x.
    """

    Background = background_current(x, p0, p1, p2, pinv) 
    peak = peak_curent(x, res, gamma, Isym, Iasym)
              
    return peak + Background

def peak_curent(x, res, gamma, Isym, Iasym):
    """
    Peak current function for ESR fitting.
    Coppied from Jose's Matlab code.

    peak_sym = peak0/(1 + (x - res)**2 / gamma**2)
    peak_asym = 2 * peak1 * (x - res) / (2 * (x - res)**2 + gamma**2)
    return peak_sym + peak_asym
    """
    
    y = 2*(x - res)/gamma 
    # don't know why the factor of 2 is needed, but it is in the original code

    nominator = Isym + 2 * Iasym * y
    denominator = y**2 + 1
    peak = nominator / denominator
    return peak

def background_current(x, p0, p1, p2, pinv):
    """
    Polynomial background current function for ESR fitting.
    Coppied from Jose's Matlab code.

    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    """
    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    return I0