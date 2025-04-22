#!/usr/bin/env python
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.append('/home/matyas/Programs/PyTimeESR')
import PyTimeESR

def Resonance(
        dyn_dict: dict,
        ham_dict: dict, 
        home_path: str,
        code_path: str, 
        step: float = .1,
        bounds: tuple = (16., 18.),
        method: str = 'brent',
        outfile: str = None, 
        code_version: str = 'bessel',
        args: dict = {}):
    """
    Run a scan of frequency and then an optimizer on the results. 
    
    Parameters:
    ----------
    dyn_dict : dict
        Dictionary containing dynamics parameters.
    ham_dict : dict
        Dictionary containing Hamiltonian parameters.
    home_path : str
        Path to the home directory where the simulation will be run.
    code_path : str
        Path to the TimeESR executable.
    step : float
        Step size for the frequency scan.
    bounds : tuple
        Frequency bounds for the scan.
    method : str
        Optimization method to be used by scipy.optimize.minimize.
    outfile : str, optional
        Path to the output file where the results will be saved.
    args : dict, optional
        Additional arguments to be passed to the scipy.optimize.minimize.

    Returns:
    -------
    minimizer : OptimizeResult
        Result of the optimization.
    frequencies : ndarray
        Array of frequencies used in the scan. 
    dc : ndarray
        Array of results from the scan.
    """    

    frequencies, amplitude = scan_freq(
        dyn_dict, ham_dict, home_path, code_path, step=step,
        bounds=bounds, outfile=outfile, code_version=code_version)
    

    # Find the index of the minimum value in the dc array
    min_index = np.argmin(amplitude)
    
    # Get the corresponding frequency
    bracket = (
        frequencies[min_index-1], frequencies[min_index], frequencies[min_index+1])
    
    minimizer = minimize_scalar(
        run_freq, bracket=bracket, 
        args=(dyn_dict, ham_dict, 
        home_path, code_path, outfile, code_version),
        method=method, **args)
    
    return minimizer, frequencies, amplitude
         
                         
def scan_freq(
        dyn_dict: dict,
        ham_dict: dict, 
        home_path: str,
        code_path: str, 
        step: float = .1,
        bounds: tuple = (16., 18.),
        outfile: str = None,
        code_version: str = 'bessel'):
    """
    Scan the resonance frequency for a given system. 
    
    Parameters:
    ----------
    dyn_dict : dict
        Dictionary containing dynamics parameters.
    ham_dict : dict
        Dictionary containing Hamiltonian parameters.
    home_path : str
        Path to the home directory where the simulation will be run.
    code_path : str
        Path to the TimeESR executable.
    step : float
        Step size for the frequency scan.
    bounds : tuple
        Frequency bounds for the scan.
    outfile : str, optional
        Path to the output file where the results will be saved.
    args : dict, optional
        Additional arguments to be passed to the scipy.optimize.minimize.
    """        
    # Create a list of frequencies to scan
    frequencies = np.arange(bounds[0], bounds[1], step)
    amplitude = np.empty_like(frequencies)

    # Run the simulation for each frequency
    for i, freq in enumerate(frequencies):
        amplitude[i] = run_freq(freq, dyn_dict, ham_dict, home_path, code_path, outfile, code_version=code_version)

    return frequencies, amplitude


def run_freq(freq: float, 
             dyn_dict: dict, 
             ham_dict: dict, 
             home_path: str, 
             code_path: str, 
             outfile: str = None, 
             code_version: str = 'bessel',
             treat_output: callable = treat_output):
    """
    Run the simulation with a given frequency and return the result.
    
    Parameters:
    ----------
    freq : float
        Frequency to be used in the simulation.
    dyn_dict : dict
        Dictionary containing dynamics parameters.
    ham_dict : dict
        Dictionary containing Hamiltonian parameters.
    home_path : str
        Path to the home directory where the simulation will be run.
    code_path : str
        Path to the TimeESR executable.
    outfile : str, optional
        Path to the output file where the results will be saved. 
    
    Returns:
    -------
    result : float
        Result of the simulation.
    """
    
    # Update the frequency in the Hamiltonian dictionary
    dyn_dict['intervals'][0]['freq'][0]['Frequency'] = freq
    
    # Create a Simulation object
    Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, home_path, 
            code_path, code_version=code_version)
    
    # Run the simulation
    Sim.run()

    _, amplitude = PyTimeESR.resonance.find_peak_freq(
        'POPULATIONS.dat', 2)
    
    if treat_output is not None:
        treat_output(Sim, outfile)
    
    return -amplitude 
    

def treat_output(Sim, outfile):
    dyn_dict = Sim.Dyn.params
    freq = dyn_dict['intervals'][0]['freq'][0]['Frequency']
    # rename POPULATIONS.dat
    #if dyn_dict['population']:
    #    os.copy('POPULATIONS.dat', f'POPULATIONS_{freq:.2f}.dat')

    if outfile is None:
        return
    
    dc = Sim.results_dict['DC']
                
    # Save the results to the specified output file
    if dyn_dict['population']:
        p, h = find_resonant_frequency(POPULATIONS.dat', 1)
        
        if not os.path.exists(outfile):
            with open(outfile, 'w') as f:
                f.write("# frequency, DC, Rabi frequency, Rabi Amplitude\n")

        with open(outfile, 'a') as f:
            f.write(f"{freq}, {dc}, {p}, {h}\n")
    else:
        if not os.path.exists(outfile):
            with open(outfile, 'w') as f:
                f.write("# frequency, DC\n")
        with open(outfile, 'a') as f:
            f.write(f"{freq}, {dc}\n")
        

if __name__ == '__main__':
    ham_dict = PyTimeESR.default_ham
    dyn_dict = PyTimeESR.default_dyn
    
    #ham_dict['eps_QD'] = -5.0
    #ham_dict['N_plot'] = 100

    dyn_dict['Ntime'] = 150000
    #dyn_dict['intervals'][0]['freq'][0]['Frequency'] = 17.02717787
    #dyn_dict['gamma0'] = [.001, .0005]
    #dyn_dict['gamma1'] = [.0, .0001]
    #dyn_dict['gammaC'] = 0.001

    #dyn_dict['biases'][0]['bias'] = [3., -3.]
    #dyn_dict['Temperature'] = 0.05
    dyn_dict['density_matrix'] = True
    dyn_dict['redimension'] = True
    dyn_dict['Nd'] = 3

    dyn_dict['population'] = True

    home_path = os.getcwd()
    code_path = '/home/matyas/Programs/TimeESR/src'
    
    Resonance(dyn_dict, ham_dict, home_path, code_path, 
        outfile='resonance.dat', bounds=(16., 18.), step = .1, 
        args={'tol': 1e-5}, code_version='standard')