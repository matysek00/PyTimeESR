import os

import numpy as np
from typing import Union

from .F90_input import F90Input
from .default_inputs import *


class Dynamics(F90Input):
    """Writer for the dynamics input. 

    Args
    -----
    dynamics_dict: dict, str
        Dictionary with the dynamics input or a file with the dynamics input.
    code_version: str
        Version of the TimeESR code to use 'standart' (default) or 'bessel'.
    line_lenght: int
        Length of a header line. Default is 80.
    padding_lenght: int
        Length of a line padding. Default is 40.
    """

    def __init__(self, dynamics_dict: Union[dict,str], code_version: str = 'standart',
                 line_lenght = 80, padding_lenght = 40):
        super(Dynamics, self).__init__(line_lenght, padding_lenght)
        

        assert code_version in ['bessel', 'standart'], \
            f"Code version {code_version} is not supported. Use 'bessel' or 'standart'."
        
        self.code_version = code_version
        self.dyn_keys = dyn_keys.copy() 

        if code_version == 'standart':
            # remove keys that are not used in the standart version
            self.dyn_keys.pop('use_bessel')
            self.dyn_keys.pop('bessel_amplitude')
            self.dyn_keys.pop('p_max')
            self.dyn_keys.pop('n_max')

        assert isinstance(dynamics_dict, (dict, str)), "Dynamics input should be a dictionary or a file name"
        if isinstance(dynamics_dict, str):
            assert os.path.exists(dynamics_dict), f"Dynamics input file {dynamics_dict} does not exist"
            self.params = self.load_input(dynamics_dict, code_version)
        else:
            self.params = dynamics_dict
        
        # check that the dictionary has the required keys with the correct types
        self.check_dictionary(self.params, self.dyn_keys, 'Dynamics input')

        N_interval = self.params['N_interval']
        Nfreq = self.params['Nfreq']
        Nbias = self.params['Nbias']

        assert N_interval > 0, "Number of intervals should be greater than 0"
        assert Nfreq > 0, "Number of frequencies should be greater than 0"
        
        assert len(self.params['intervals']) == N_interval, "Number of intervals should match the number of intervals in the dictionary"
        assert len(self.params['biases']) == Nbias, "Number of biases should match the number of biases in the dictionary"

        for i in range(N_interval): 
            # check that the intervals have the required keys with the correct types
            self.check_dictionary(self.params['intervals'][i], dyn_keys_interval, 'Dynamics input - Intervals')

            assert len(self.params['intervals'][i]['freq']) == Nfreq, "Number of frequencies should match the number of frequencies in the dictionary"
            
            for j in range(Nfreq):
                # check that the frequencies have the required keys with the correct types
                self.check_dictionary(self.params['intervals'][i]['freq'][j], dyn_keys_freq, 'Dynamics input - Frequencies')

        for i in range(Nbias):
            # check that the biases have the required keys with the correct types
            self.check_dictionary(self.params['biases'][i], dyn_keys_bias, 'Dynamics input - Biases')

        intervals = self.params['intervals']
        assert all(intervals[i]['t0'] < intervals[i]['tf'] for i in range(N_interval)), "Interval times should be in increasing order"
        assert all(intervals[i]['t0'] == intervals[i-1]['tf'] for i in range(1, N_interval)), "Interval times should be continuous"
        assert intervals[0]['t0'] == 0, "First interval should start at 0"
        assert intervals[-1]['tf'] == self.params['t_final'], "Last interval should end at t_final"
        
        if code_version == 'standart':
            return
        
        if  self.params['use_bessel']:
            assert Nfreq == 1, "Number of frequencies should be 1 for Bessel function"
            assert N_interval == 1, "Number of intervals should be 1 for Bessel function"
            assert self.params['p_max'] > self.params['n_max'], "Max order of Bessel function should be greater than max frequency"
            assert self.params['n_max'] >= 0, "Max frequency of Bessel function should be greater than or equal to 0"


    def write_input(self,):
        """Write the dynamics into a string.

        Returns
        -------
        str: Dynamics input string.
        """

        input_string = self.create_header('', '*')
        input_string += self.create_header('TimeESR Input', '*')
        input_string += self.create_header('', '*')
        input_string += self.input_line(self.params['Ntime'], 'Number of time points')
        input_string += self.input_line([self.params['t_initial'], self.params['t_final']], 'Initial and final time (ns)')
        
        N_interval = self.params['N_interval']
        Nfreq = self.params['Nfreq']
        input_string += self.create_header('Pulse definition block', '-')
        input_string += self.input_line(N_interval, 'Number of intervals')
        input_string += self.input_line(Nfreq, 'Number of frequencies per interval')
        for i in range(N_interval):
            interval = self.params['intervals'][i]
            input_string += self.input_line([interval['t0'],interval['tf']], f'{i} - times for pulses (ns)')

            for j in range(Nfreq):
                freq = interval['freq'][j]
                input_string += self.input_line(freq['Amplitude'], f'{i} - F{j} amplitude')
                input_string += self.input_line(freq['Frequency'], f'{i} - F{j} frequency (GHz)')
            input_string += self.input_line(interval['Phase'], f'{i} - phase (radians)')
        
        input_string += self.create_header('Electrode set-up ', '-')
        input_string += self.input_line(self.params['gamma0'][0], 'gamma_R_0= 2*pi*W_R_0*W_R_0*rho (meV)')
        input_string += self.input_line(self.params['gamma0'][1], 'gamma_L_0 2*pi*W_L_0*W_L_0*rho (meV)')
        input_string += self.input_line(self.params['gamma1'][0], 'gamma_R_1= 2*pi*W_R_0*W_R_1*rho (meV)')
        input_string += self.input_line(self.params['gamma1'][1], 'gamma_R_1= 2*pi*W_L_0*W_L_1*rho (meV)')
        input_string += self.input_line(self.params['cutoff'], 'Cutoff for integral Lambshift (meV)')
        input_string += self.input_line(self.params['gammaC'], 'Broadening of Green\'s function (meV)')
        input_string += self.input_line(self.params['integral_points'], 'Number of points for  I11 and I21')

        Nbias = self.params['Nbias']
        input_string += self.create_header('Bias, temperature, and spin polarization', '-')
        input_string += self.input_line(Nbias, 'Number of biases')
        for i in range(Nbias):
            bias = self.params['biases'][i]
            input_string += self.input_line(bias['bias'][0], f'{i} - Right electrode bias (mV)')
            input_string += self.input_line(bias['bias'][1], f'{i} - Left electrode bias (mV)')
            input_string += self.input_line(bias['b_time'], f'{i} - Duration of the pulse (ns)')
        
        input_string += self.input_line(self.params['Temperature'], 'Temperature (K)')
        input_string += self.input_line(self.params['spin_polarization'][0], 'Right electrode spin polarization')
        input_string += self.input_line(self.params['spin_polarization'][1], 'Left electrode spin polarization')
        input_string += self.input_line(self.params['Electrode'], 'Current measurement: 0 is left and 1 is right electrode')

        if self.code_version == 'bessel':
            input_string += self.create_header('Bessel function', '-')
            input_string += self.input_line(self.params['use_bessel'], 'Use Bessel function')
            input_string += self.input_line(self.params['bessel_amplitude'][0], 'B_R strengt of the time depenndet pulse for right electrode')
            input_string += self.input_line(self.params['bessel_amplitude'][1], 'B_L strengt of the time depenndet pulse for left electrode')
            input_string += self.input_line(self.params['p_max'], 'Max order of Bessel function in both directions')
            input_string += self.input_line(self.params['n_max'], 'Max frequency of Bessel function in both directions')
        
        input_string += self.create_header('Output', '-')
        input_string += self.input_line(self.params['population'], 'write POPULATION.dat')
        input_string += self.input_line(self.params['density_matrix'], 'write rho.dat')
        input_string += self.input_line(self.params['output_file'], 'file name for the time-dependent current')
        input_string += self.input_line(self.params['output_fourier'], 'file name for the Fourier transform of the current')
        input_string += self.input_line(self.params['output_ESR'], 'file for DC current')

        input_string += self.create_header('read previous run (DANGEROUS)', '-')
        input_string += self.input_line(self.params['runs'], 'Keep to false to avoid reading previous runs')

        input_string += self.create_header('other options', '-')
        input_string += self.input_line(self.params['spindyn'], 'Plot spin dynamics')
        input_string += self.input_line(self.params['redimension'], 'Redimension the Hamiltonian')
        input_string += self.input_line(self.params['Nd'], 'Redimension number - must include empty states for non-zero current')\
        
        input_string += self.create_header('', '*')
        input_string += self.create_header('End of input', '*')
        input_string += self.create_header('', '*')
        return input_string
    
    @staticmethod
    def load_input(input_file: str, code_version: str = 'standart'):
        """Load the dynamics from a file.
        
        Args
        -----
        input_file: str
            Dynamics input file
        code_version: str
            Version of the TimeESR code to use standart (default) or bessel.

        Returns
        -------
        params: dict
            Dictionary with the dynamics input
        """
        print('WARNING: Loading dynamics from file is not tested')

        infile = open(input_file, 'r')
        params = {}

        _ = infile.readline()
        _ = infile.readline()
        _ = infile.readline()

        params['Ntime'] = int(infile.readline().split()[0])
        
        times = [float(x) for x in infile.readline().split()[:2]]
        params['t_initial'] = float(times[0])
        params['t_final'] = float(times[1])
        
        _ = infile.readline()
        params['N_interval'] = int(infile.readline().split()[0])
        params['Nfreq'] = int(infile.readline().split()[0])

        params['intervals'] = []
        for _ in range(params['N_interval']):
            interval = {}
            
            times = [float(x) for x in infile.readline().split()[:2]]
            interval['t0'] = float(times[0])
            interval['tf'] = float(times[1])
            interval['freq'] = []
            for _ in range(params['Nfreq']):
                freq = {}
                freq['Amplitude'] = float(infile.readline().split()[0])
                freq['Frequency'] = float(infile.readline().split()[0])
                interval['freq'].append(freq)
            interval['Phase'] = float(infile.readline().split()[0])
            params['intervals'].append(interval)

        _ = infile.readline()
        params['gamma0'] = [float(infile.readline().split()[0]), 
                            float(infile.readline().split()[0])]
        params['gamma1'] = [float(infile.readline().split()[0]), 
                            float(infile.readline().split()[0])]
        params['cutoff'] = float(infile.readline().split()[0])
        params['gammaC'] = float(infile.readline().split()[0])
        params['integral_points'] = int(infile.readline().split()[0])
        
        params['biases'] = []
        _ = infile.readline()
        params['Nbias'] = int(infile.readline().split()[0])
        for _ in range(params['Nbias']):
            bias = {}
            bias['bias'] = [float(infile.readline().split()[0]), 
                            float(infile.readline().split()[0])]
            bias['b_time'] = float(infile.readline().split()[0])
            params['biases'].append(bias)

        params['Temperature'] = float(infile.readline().split()[0])
        params['spin_polarization'] = [float(infile.readline().split()[0]), 
                                       float(infile.readline().split()[0])]
        params['Electrode'] = int(infile.readline().split()[0])

        if code_version == 'bessel':
            _ = infile.readline()
            params['use_bessel'] = F90Input.string2bool(infile.readline().split()[0])
            params['bessel_amplitude'] = [float(infile.readline().split()[0]), 
                                          float(infile.readline().split()[0])]
            params['p_max'] = int(infile.readline().split()[0])
            params['n_max'] = int(infile.readline().split()[0])

        _ = infile.readline()
        params['population'] = F90Input.string2bool(infile.readline().split()[0])
        params['density_matrix'] = F90Input.string2bool(infile.readline().split()[0])
        params['output_file'] = infile.readline().split()[0]
        params['output_fourier'] = infile.readline().split()[0]
        params['output_ESR'] = infile.readline().split()[0]

        _ = infile.readline()
        params['runs'] = F90Input.string2bool(infile.readline().split()[0])

        _ = infile.readline()
        params['spindyn'] = F90Input.string2bool(infile.readline().split()[0])
        params['redimension'] = F90Input.string2bool(infile.readline().split()[0])
        params['Nd'] = int(infile.readline().split()[0])
        infile.close()
        
        return params
        
    def create_output_dict(self):
        
        output_dict = {}
        output_dict['current'] = self.params['output_file']
        output_dict['fourier'] = self.params['output_fourier']
        output_dict['DC'] = self.params['output_ESR']

        if self.params['population']:
            output_dict['population'] = 'POPULATIONS.dat'
        if self.params['density_matrix']:
            output_dict['density_matrix'] = 'RHO.dat'
        if self.params['spindyn']:
            output_dict['spin_dynamics'] = 'SpinDynamics.dat'

        return output_dict
    
    def load_output(self, output_dict):
        result_dict = {}
        self.results_dict['DC'] = np.loadtxt(self.output_dict['DC'])
