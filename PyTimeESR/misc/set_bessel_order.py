import numpy  as np

from scipy.special import jv

def find_bessel_order(Vrf: float, omega: float, pmax: int = 1e6, tresh=1e-4): 
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
    Vrf *= 1e-3 # eV
    omega *= 1e9 # Hz 

    z = Vrf/omega
    
    ps = np.arange(0,pmax, 1, dtype=int)
    J = np.abs(jv(ps,z))
    above = J < 1e-5

    for p in ps[::-1]:
        if not above[p]:
            pbes = p+1
            break
    if pbes >= pmax:
        Warning(f'Provided pmax={pmax:d} is insuficient.')
    
    return pbes
    