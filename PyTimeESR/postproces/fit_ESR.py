import numpy as np
from scipy.optimize import curve_fit



def create_fit_esr(x, y, Zeeman, lb=None, ub=None, maxfev=None, size = 1.,tol=1e-4, weight_scale=0):
    """Probably redundant funcition
   
    p0 = (res, gamma, Isym, Iasym, p0, p1, p2, pinv)
    
    Ips = Isym/(1 + 4 * (x - res)**2 / gamma**2)
    Ipa = 2 * Iasym * (x - res) / gamma / (1 + 4 * (x - res)**2/ gamma**2)
    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    I = I0 + Ips + Ipa
    weighscale: if >0, apply a weight that penalizes points away from guessed resonance
        to put more emphasis on fitting the peak region.
    """
    ### NEED To FIT ALL Peaks at SIMULTANEOUSLY
    p0 = guess_p0(x, y, Zeeman, size)
    use = np.abs(x - Zeeman) < size
    
    var = np.max(y[use]) - np.min(y[use])
    
    if var < tol:
        opt = np.zeros_like(p0)
        opt[0] = Zeeman
        return np.zeros_like(p0), np.zeros_like(p0)
    try:
        popt, pcov = curve_fit(ESR_fit_fun, x[use], y[use], p0=p0, bounds=(lb, ub), maxfev=maxfev,)       
    except RuntimeError:
        opt = np.zeros_like(p0)
        opt[0] = Zeeman
        return np.zeros_like(p0), np.zeros_like(p0)

    return popt, np.sqrt(np.diag(pcov))

def ESR_fit_fun(x, res, gamma, Isym, Iasym, p0, p1, p2, pinv):
    """
    ESR fitting function.
    
    Ips = Isym/(1 + (x - res)**2 / gamma**2)
    Ipa = 2 * Iasym * (x - res) / (2 * (x - res)**2 + gamma**2)
    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    I = I0 + Ips + Ipa

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

    nominator = Isym + Iasym * y
    denominator = y**2 + 1
    peak = nominator / denominator
    return peak

def background_current(x, p0, p1, p2, pinv,):
    """
    Polynomial background current function for ESR fitting.
    Coppied from Jose's Matlab code.

    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    """
    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    return I0

def guess_p0(x, y, x0, size):
    
    b = (y[-1]-y[0])/(x[-1]-x[0])
    a = y[0] - b*x[0]
    
    use = np.abs(x - x0) < size

    y = y[use]
    x = x[use]

    y -= b*x + a

    Is = y[np.argmin(np.abs(x-x0))]


    argmin = np.argmin(y)
    argmax = np.argmax(y)
    
    gamma = x[argmax] - x[argmin] # assuming Is = 0, but should be close
    Ia = y[argmax] - y[argmin]
    
    Ia = gamma/np.abs(gamma)*Ia # fix sign
    gamma = np.abs(gamma)


    p0 = [x0, gamma, Is, Ia, a, b , 0., 0.]
    upper_bound = [x0 + size/2, 1.,  1.1*np.abs(Ia),  1.1*np.abs(Ia), np.inf, np.inf, np.inf, np.inf]
    lower_bound = [x0 - size/2, 0., -1.1*np.abs(Ia), -1.1*np.abs(Ia), np.inf, np.inf, np.inf, np.inf]
    return p0