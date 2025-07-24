import numpy as np 

def w0c(x: float, A: float, phi: float) -> float:
    acos = np.arccos(x)
    sq = np.sqrt(1-x**2)

    w0 = 1/np.pi*(
        (1+A**2/2)*acos
        + 2*A*sq*np.cos(phi) 
        + A**2/2*x*sq*np.cos(2*phi)
    )
    return w0

def w1c(x: float, A: float, phi: float) -> float:
    acos = np.arccos(x)
    sq = np.sqrt(1-x**2)
    phase = np.exp(phi*1.j)

    w1  = 1/np.pi*(
          A*phase*acos 
        + A*phase*x*sq 
        + sq*( 1 + (A**2)/3 * (2*np.cos(phi)**2 + np.sin(phi)**2))
        + (A**2)/3*sq*(x**2)*(phase**2))
    return w1 

def w0s(x: float, A: float, phi: float) -> float:
    w0 = (1+A**2/2) - w0c(x,A,phi)
    return w0

def w1s(x: float, A: float, phi: float) -> float:
    w0 = A - w1c(x,A,phi)
    return w0

t_four_comps_cot = {
    0: w0c,
    1: w1c
}

t_four_comps_seq = {
    0: w0c,
    1: w1c
}

def four_background_current(
        I0: float, A: float, x: float, phi:float, 
        m: int=0, regime: str='seq') -> float:
    """Estimate backround current (I_back).

    Args
    -----------
    I0: float
        Current for A=0, B=0
    A: float
        Tunel driving amplitude 
    x: float
        (epsilon_QD - Ef)/B
    phi: float
        phase betwee A and B driving 
    m: int
        fourier component of the current (default 0)
    regime: str
        'seq' (default), or 'cot'

    Returns:
        I: float
    """

    if regime not in ['seq', 'cot']: 
        raise ValueError('regime must be \'cot\' or \'seq\'')
    t_four_comps = t_four_comps_seq if regime == 'seq' else t_four_comps_cot
    
    x = np.array([x]) if type(x) is not np.array else x 
    x[x>1.] = 1. 

    tm = t_four_comps[m](x,A,phi)
    
    return I0*tm


def four_peak_current(domega: float, A: float, x: float, 
        Gamma0: float, phi: float, n: int=1,
        m: int=0, regime: str='seq') -> float:
    """Estimate current Peak current (Ip).

    Args
    -----------
    domega: float
        n*omega - Zeman energy
    A: float
        Tunel driving amplitude 
    x: float
        (epsilon_QD - Ef)/B
    Gamma0: float
        Re(\sum_{v=2} \Gamma_{1vv1}(A=0, B=0))
    phi: float
        phase betwee A and B driving 
    n: int
        resonance number default 1
    m: int
        fourier component of the current (default 0)
    regime: str
        'seq' (default), or 'cot'

    Returns:
        I_m: float
    """
    
    if regime not in ['seq', 'cot']: 
        raise ValueError('regime must be \'cot\' or \'seq\'')
    t_four_comps = t_four_comps_seq if regime == 'seq' else t_four_comps_cot
    
    x = np.array([x]) if type(x) is not np.array else x 
    x[x>1.] = 1. 

    t0 = t_four_comps[0](x,A,phi)
    tn = t_four_comps[n](x,A,phi)
    tm = t_four_comps[m](x,A,phi)
    
    I = 2*Gamma0 * tn*tm*t0 / (t0**2 + (domega/Gamma0)**2)
    return I