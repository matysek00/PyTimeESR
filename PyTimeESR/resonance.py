import os
import numpy as np

from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar

from .pytimeesr import Simulation

def find_resonant_frequency(fil, col, freq_bound = (1e-6, 1.)):
    """Find dominant freqeuncy and its amplitude in a POPULATIONS.dat file.

    Parameters
    ----------
    fil : str
        Path to the POPULATIONS.dat file.
    col : int
        Column number to analyze.
    freq_bound : tuple, optional
        Frequency bounds for the analysis. Default is (1e-6, 1.).
    
    Returns
    -------
    tuple
        Dominant frequency and its amplitude.
    """

    # Load the data from the file
    pop = np.loadtxt(fil, usecols=col)
    time = np.loadtxt(fil, usecols=0)
    
    # Fourier transform the data
    dt = (time[-1] - time[0])/len(time)
    freqs = fftfreq(len(pop), dt)

    # Take the absolute value of the real component of the FFT
    # in case decay would be more pronounced than oscillation
    # and remove the zero frequency component
    fpop = np.abs(np.real(fft(pop- np.mean(pop))))
    
    # Chose frequency range
    frange_plus = (freqs < freq_bound[1])*(freqs > freq_bound[0])*(freqs >= 0)
    frange_mins = (freqs >-freq_bound[1])*(freqs <-freq_bound[0])*(freqs < 0)
    fpopt = np.concatenate((fpop[frange_mins],fpop[frange_plus]))
    freqst = np.concatenate((freqs[frange_mins], freqs[frange_plus]))
 
    # Find the peaks in the Fourier transform
    peaks = find_peaks(fpopt)[0]
    if len(peaks) == 0:
        # No peaks found, return zero frequency and amplitude
        return 0, 0. 
    
    # Find the peak with the highest amplitude
    heights = fpopt[peaks]
    pm = peaks[np.argmax(heights)]

    return freqst[pm], np.max(heights)


def treat_output(Sim,
                  outfile: str = None,
                  ):
    """
    Treat the output of the simulation and save it to a file.
    
    Parameters:
    ----------
    Sim : Simulation
        Simulation object containing the results.
    outfile : str
        Path to the output file where the results will be saved.
    """
    dyn_dict = Sim.Dyn.params
    freq = dyn_dict['intervals'][0]['freq'][0]['Frequency']

    if outfile is None:
        return
    
    dc = Sim.results_dict['DC']
                
    # Save the results to the specified output file
    if dyn_dict['population']:
        rab_freq, rab_amp = find_resonant_frequency('POPULATIONS.dat', 1)
        
        if not os.path.exists(outfile):
            with open(outfile, 'w') as f:
                f.write(
                 "# frequency, DC, Rabi frequency, Rabi Amplitude\n")

        with open(outfile, 'a') as f:
            f.write(f"{freq}, {dc}, {rab_freq}, {rab_amp}\n")
    else:
        if not os.path.exists(outfile):
            with open(outfile, 'w') as f:
                f.write("# frequency, DC\n")
        with open(outfile, 'a') as f:
            f.write(f"{freq}, {dc}\n")


def run_freq(freq: float, 
             dyn_dict: dict, 
             ham_dict: dict, 
             home_path: str, 
             code_path: str, 
             outfile: str = None, 
             code_version: str = 'standart',
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
    Sim = Simulation(ham_dict, dyn_dict, home_path, 
            code_path, code_version=code_version)
    
    # Run the simulation
    Sim.run()

    _, amplitude = find_resonant_frequency(
        'POPULATIONS.dat', 2)
    
    if treat_output is not None:
        treat_output(Sim, outfile)
    
    return -amplitude 


def scan_freq(
        dyn_dict: dict,
        ham_dict: dict, 
        home_path: str,
        code_path: str, 
        step: float = .1,
        bounds: tuple = (16., 18.),
        outfile: str = None,
        code_version: str = 'standart'):
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
    code_version : str
        Version of the TimeESR code to use 'standart' (default) or 'bessel'.
    """        
    # Create a list of frequencies to scan
    frequencies = np.arange(bounds[0], bounds[1], step)
    amplitude = np.empty_like(frequencies)

    # Run the simulation for each frequency
    for i, freq in enumerate(frequencies):
        amplitude[i] = run_freq(freq, dyn_dict, ham_dict, home_path, code_path, outfile, code_version=code_version)

    return frequencies, amplitude


def Resonance(
        dyn_dict: dict,
        ham_dict: dict, 
        home_path: str,
        code_path: str, 
        step: float = .1,
        bounds: tuple = (16., 18.),
        method: str = 'brent',
        outfile: str = None, 
        code_version: str = 'standart',
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
    code_version : str
        Version of the TimeESR code to use 'standart' (default) or 'bessel'.
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
    
    # check if the minimum is at the edge of the array
    if min_index == 0 or min_index == len(amplitude)-1:
        print("WARNING: Minimum is at the edge of the boundary")
        return None, frequencies, amplitude
    
    # Check if amplitude is constant
    if np.all(amplitude == amplitude[0]):
        print("WARNING: Amplitude is constant accross the range")
        return None, frequencies, amplitude    

    # Check if the bracket is valid
    bracket_works = ((
        amplitude[min_index] < amplitude[min_index-1]) and(
        amplitude[min_index] < amplitude[min_index+1]))

    # Get the corresponding frequency
    bracket = (
        frequencies[min_index-1], frequencies[min_index], frequencies[min_index+1])

    if not bracket_works:
        print("WARNING: No minimum found")
        return None, frequencies, amplitude

    minimizer = minimize_scalar(
        run_freq, bracket=bracket, 
        args=(dyn_dict, ham_dict, 
        home_path, code_path, outfile, code_version),
        method=method, **args)
    
    return minimizer, frequencies, amplitude