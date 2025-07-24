from scipy.optimize import minimize


def Iset_step(gL, Itarget, Sim):

    Sim.Dyn.params['gamma0'][1] = gL
    Sim.run()
    Sim.load_output()
    DC = Sim.results_dict['DC']
    print(DC)
    return (DC-Itarget)**2


def Iset(Sim, Itarget):
    # initial guess
    gL = Sim.Dyn.params['gamma0'][1]
    
    # turn off driving 
    Sim.Dyn.params['A'] = [0+0j,0+0j]
    Sim.Dyn.params['bessel_amplitude'] = [0,0]
    Sim.Dyn.params['n_max'] = 3
    Sim.Dyn.params['p_max'] = 5

    res = minimize(Iset_step, gL, (Itarget, Sim))
    return res
