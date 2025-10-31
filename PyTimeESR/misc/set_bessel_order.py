import numpy  as np

from scipy.special import jv

def find_bessel_order(Vrf: float, omega: float, pmax: int = 1e6, tresh=.001): 
    """

    Parameters:
    -----------
        Vrf (float):
            in meV 
        omega (float):
            in GHz
        pmax (int):
            maximum possible
    """
    hbar = 6.582119569e-16
    Vrf *= 1e-3 # eV
    omega *= 1e9 # Hz 

    z = Vrf/(2*np.pi*hbar*omega)
    
    ps = np.arange(0,pmax, 1, dtype=int)
    J = np.abs(jv(ps,z))
    
    csJ = np.cumsum(J)
    csJ /= csJ[-1]
    above = csJ > 1-tresh
    
    pbes = pmax + 1 
    for p in ps[::]:
    
        if above[p]:
            pbes = p+1
            break 
    
    if pbes >= pmax:
        Warning(f'Provided pmax={pmax:d} is insuficient.')
    
    pbes += 1 # to include shift from Adriving 
    return pbes
    