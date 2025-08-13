from scipy.optimize import root_scalar
import numpy as np

def Iset(Sim, 
         Itarget: float, 
         outfile: str = None, 
         frequency: float = None, 
         init_step: float = None,
         method: str = 'brentq',
         ratio: float = 1.):
    """
    Set the current in the simulation by adjusting the damping parameter gamma0.
    Bound by 0, and the right damping parameter gamma0[0].
    Similiar (with different optimization method) to the Iset function in the Fortran code.
    
    Parameters
    ----------
    Sim : Simulation object
        The simulation object containing the dynamic parameters.
    Itarget : float
        The target current to set in the simulation.
    outfile : str, optional
        Where to dump Fortran output, by default None. 
    frequency : float, optional
        The frequency of the driving force, by default None (no driving).
        WARNING: A and B drive parameters are set to zero if frequency is None.
    init_step : float, optional
        Initial step size for the root finding method, by default None.
    method : str, optional
        The method to use for root finding, by default 'brentq'.
    ratio : float, optional
        Ratio to scale the upper bound of gL/gR, by default 1.0.
    
    Returns
    -------
    res : OptimizeResult
        If suceesful, the result of the optimization containing the root and convergence status 
        If not successful, a dictionary with the root and a message indicating failure.
    """
    
    # initial guess
    gL = Sim.Dyn.params['gamma0'][1]
    gR = Sim.Dyn.params['gamma0'][0]
    bounds = (1e-8, ratio*gR)
    x1 = None if init_step is None else gL+init_step # secondary guess
    
    Itarget = np.abs(Itarget)  # ensure target current is positive
    
    if frequency is None:
        # turn off driving 
        print('WARNING: No frequency provided, setting driving parameters to zero.', flush=True)
        Sim.Dyn.params['A'] = [0+0j,0+0j]
        Sim.Dyn.params['bessel_amplitude'] = [0,0]
        Sim.Dyn.params['n_max'] = 3
        Sim.Dyn.params['p_max'] = 5
    else: 
        Sim.Dyn.params['frequency'] = frequency

    #res = minimize(Iset_step, gL, (Itarget, Sim, outfile), bounds = bounds, options=options)
    try:
        res = root_scalar(Iset_step, args = (Itarget, Sim, outfile), 
            method=method, bracket=bounds, x0=gL, x1=x1, xtol=1e-8)
    except: 
        res = {'converged': False, 'root': gL, 'message': 'Failed with method: ' + method}
    return res


def Iset_step(gL, Itarget, Sim, outfile=None):

    Sim.Dyn.params['gamma0'][1] = gL
    Sim.run(outfile)
    Sim.load_output()
    DC = np.abs(Sim.results_dict['DC'])

    print('Iset_step: gammaL =', gL, 'DC =', DC, 'Itarget =', Itarget, flush=True)
    return DC-Itarget

def Feedback(Sim, 
         Itarget: float, 
         outfile: str = None, 
         tol: float = 1e-6,
         n_steps: int = 1e3):

    gl  = Sim.Dyn.params['gamma0'][1]
    gr  = Sim.Dyn.params['gamma0'][0]
    
    Itarget = np.abs(Itarget)  # ensure target current is positive

    for i in range(int(n_steps)):
        # usign Iset=0 for Iset_step funciont
        I = Iset_step(gl, 0.0, Sim, outfile) 
        dI = I - Itarget
        
        print(f'Step {i+1}: Current = {I}, Target = {Itarget}, dI = {dI}, gamma0 = {gl}')

        if np.abs(dI) < tol:
            print('Feedback converged: Current is within tolerance of target.')
            return {'converged': True, 'root': gl, 'message': 'Feedback converged'}
        
        gl += (5 * np.sign(dI) * tol * (dI - tol/10)/(I+Itarget) 
            * (1+(-1)**(np.int(1 - .1*np.sign(dI)))))
        gl += (5 * np.sign(dI) * tol * (-dI - tol/10)/(I+Itarget) 
            * (1-(-1)**(np.int(1 - .1*np.sign(dI)))))
        
        if (gl <= 0 or gr <= gl):
            print('Feedback failed: gamma0 bounds are not positive:', gl, gr)
            return {'converged': False, 'root': gl, 'message': 'gamma0 out of bounds'}
        
    # if we reach here, feedback did not converge
    print('Feedback did not converge within the maximum number of steps.')
    return {'converged': False, 'root': gl, 'message': 'Feedback did not converge'}