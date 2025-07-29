from scipy.optimize import minimize, root_scalar
import numpy as np

def Iset_step(gL, Itarget, Sim, outfile=None):

    Sim.Dyn.params['gamma0'][1] = gL
    Sim.run(outfile)
    Sim.load_output()
    DC = Sim.results_dict['DC']
    return DC-Itarget


def Iset(Sim, Itarget, outfile=None, frequency=None, init_step=None,
         method='brentq'):#options={}):
    # initial guess
    
    gL = Sim.Dyn.params['gamma0'][1]
    gR = Sim.Dyn.params['gamma0'][0]
    bounds = (1e-8, gR)
    x1 = None if init_step is None else gL+init_step # secondary guess
    
    if frequency is None:
        # turn off driving 
        Sim.Dyn.params['A'] = [0+0j,0+0j]
        Sim.Dyn.params['bessel_amplitude'] = [0,0]
        Sim.Dyn.params['n_max'] = 3
        Sim.Dyn.params['p_max'] = 5
    else: 
        Sim.Dyn.params['frequency'] = frequency

    #res = minimize(Iset_step, gL, (Itarget, Sim, outfile), bounds = bounds, options=options)
    res = root_scalar(Iset_step, args = (Itarget, Sim, outfile), 
        method=method, bracket=bounds, x0=gL, x1=x1, xtol=1e-8)
    return res
