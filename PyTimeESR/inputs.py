
import os
import re

from typing import Union

from .default_inputs import *


class _complex_format(str):
    """Custom string class to format complex numbers.
    """

    def format(self):
        return '({:.8f},{:.8f})'.format(self.real, self.imag)              


class F90Input(): 
    """Parent class to write fortran input.
    """
    
    fmt = {
        tfloat: '{:.8f} ',
        tcomplex: _complex_format,
        tint: '{:d} ',
        tstr: '{} ',
        tbool: '{} ',
    }
    
    def __init__(self, line_lenght, padding_lenght):
        self.line_lenght = line_lenght
        self.padding_lenght = padding_lenght

    @staticmethod
    def check_dictionary(dictionary, required_keys, dict_label):
        """Check if the dictionary has the required keys and types.
        Args
        -----
        dictionary: dict
            Dictionary to check.
        required_keys: dict
            Dictionary with the required keys and types. if the type is a list, the first element
            is the type of the list, the second element type of elememts,
            and the third element is the number of elements in the list.
        dict_label: str
            Label of the dictionary to check. to be used in the error message.
        """
        assert isinstance(dictionary, dict), "Hamiltonian input should be a dictionary"
        for key, value in required_keys.items(): 
            assert key in dictionary, f"Missing key {key} in {dict_label}"

            if not isinstance(value, tlist):
                # if not a list check if the type is correct
                assert isinstance(dictionary[key], value), f"Key {key} should be of type {value} in {dict_label}"
                continue
            
            # if itterable check if the type is correct
            assert isinstance(dictionary[key], value[0]), f"{key} should be of type {value} in {dict_label}"
            # check if the elements of the list are of the correct type
            assert all(isinstance(i, value[1]) for i in dictionary[key]), f"Elements of {key} should be of type {value} in {dict_label}"
            if value[2] is None:
                continue
            # check if the number of elements is correct if it is not fixed
            assert len(dictionary[key]) == value[2], f"Key {key} should have {value[2]} elements in {dict_label}"

    def create_header(self, header: str, seperator: str):
        """Create the header for the input file.
        Args
        -----
        header: str
            Header to add to the input file.
        seperator: str
            Seperator to add to the input file.
        lenghth: int
            Length of the header. Default is 80.
        """

        l = int((self.line_lenght - len(header))/2)
        line = seperator * l + header + seperator * l + '\n'
        return line
    
    @staticmethod
    def bool2string(value: bool):
        string = '.true.' if value else '.false.'
        return string
    
    @staticmethod
    def string2bool(value: str):
        if value == '.true.':
            return True
        elif value == '.false.':
            return False
        else:
            raise ValueError(f"Invalid boolean string: {value}. Use '.true.' or '.false.'.")
        
    @staticmethod
    def load_complex(line: str) -> complex:
        """Load a complex number from a line of the input file.
        format: (Re, Im)

        Args
        -----
        line: str
            Line from the input file containing a complex number.

        Returns
        -------
        complex: 
            Complex number loaded from the line.
        """
        string = re.split(r'[()\s,]', line)[1:3]
        z = complex(*tuple(map(float, string)))
        return z

    def input_line(self, value, comment: str = ''):
        """Create a line for the input file.
        Args
        -----
        value: str
            Value to add to the input file.
        comment: str
            Comment to add to the input file.
        """
        # turn everything into a list
        if not isinstance(value,tlist):
            value = [value]
        
        
        # convert bools to strings
        value = [self.bool2string(val) if isinstance(val, tbool) else val for val in value]
        t = type(value[0])
    
        # convert t the specified type into a generic type ttype
        ttype = tfloat if (t in tfloat) else tint if (t in tint) else tstr if (t in tstr) else tcomplex if (t in tcomplex) else None
        assert ttype is not None, f"Type {t} not supported in input line"
    
        line = ' '.join([self.fmt[ttype].format(i) for i in value])
        pad = self.padding_lenght - len(line)
        pad = 0 if pad < 0 else pad
        line += ' ' * pad + '! ' + comment + '\n'
        return line
    
    def load_output(self, output_dict: str):
        """Load the output from a file.
        
        Args
        -----
        output_dict: str
            Output file to load.
        
        Returns
        -------
        dict: Dictionary with the output.
        """
        
        return {}
    

class Hamiltonian(F90Input):
    """Write input for the Hamiltonian. 

    Args
    -----
    hamiltonian_dict: (dict, str)
        Dictionary with the Hamiltonian input or a file with the Hamiltonian input.
    """
    
    def __init__(self, hamiltonian_dict: Union[dict, str],
                  line_lenght = 80, padding_lenght = 40):
        super(Hamiltonian, self).__init__(line_lenght, padding_lenght)

        assert isinstance(hamiltonian_dict, (dict, str)),  "Hamiltonian input should be a dictionary or a file name"
        if isinstance(hamiltonian_dict, str):
            assert os.path.exists(hamiltonian_dict), f"Hamiltonian input file {hamiltonian_dict} does not exist"
            self.load_input(hamiltonian_dict)
        else:
            self.params = hamiltonian_dict
        self.check_dictionary(self.params, ham_keys, 'Hamiltonian input')
        
        Nm = self.params['Nm']
        Npairs = self.params['Npairs']
        
        assert  Nm > 0, "Number of spins should be greater than 0"
        assert Npairs >= 0, "Number of pairs should be greater than or equal to 0"
        
        assert len(self.params['Spins']) == Nm, "Number of spins should match the number of spins in the dictionary"
        assert len(self.params['pairs']) == Npairs, "Number of pairs should match the number of pairs in the dictionary"

        for i in range(Nm): 
            self.check_dictionary(self.params['Spins'][i], ham_keys_spins, 'Hamiltonian input - Spins')
        for i in range(Npairs): 
            self.check_dictionary(self.params['pairs'][i], ham_keys_pairs, 'Hamiltonian input - Pairs')
        
        S_trans = self.params['Spins'][0]['S']
        assert S_trans == .5, 'Transport electron should be a spin 1/2'
        

    def write_input(self,):
        """Write the Hamiltonian into a string.

        Returns
        -------
        str: Hamiltonian input string.
        """
        
        Nm = self.params['Nm']
        
        input_string = self.create_header('', '*')
        input_string += self.create_header('Hamiltonian Input', '*')
        input_string += self.create_header('', '*')
        input_string += self.input_line(Nm, 'Number of spins')

        for i in range(Nm):
            Spin = self.params['Spins'][i]
            input_string += self.create_header(f'Spin Properties: Spin {i}', '-')
            input_string += self.input_line(Spin['S'], 'Spin')
            input_string += self.input_line(Spin['Stephen'], 'Stephen Coefficients: B20 B22 B40 B42')
            input_string += self.input_line(Spin['Stephen_axis'], 'Stephen Axis')
            input_string += self.input_line(Spin['H'], 'Local Magnetic Field: Hx Hy Hz')
            input_string += self.input_line(Spin['G'], 'Local gyrometric factor: Gxx Gyy Gzz')
        

        Npairs = self.params['Npairs']
        input_string += self.create_header('Exchange Interactions', '-')
        input_string += self.input_line(Npairs, 'Number of connected pairs')
        for i in range(Npairs):
            pair = self.params['pairs'][i]
            input_string += self.input_line(pair['pair'], 'Connected pair: i j')
            input_string += self.input_line(pair['J_exc'], 'Exchange interaction: Jxx Jxy Jxz (GHz)')
        
        input_string += self.create_header('Electronic Interactions', '-')
        input_string += self.input_line(self.params['eps_QD'], 'QD energy (meV)')
        input_string += self.input_line(self.params['Hubbard'], 'Hubbard interaction U(meV)')

        input_string += self.create_header('Output', '-')
        input_string += self.input_line(self.params['output_file'], 'Output file for the Hamiltonian')
        input_string += self.input_line(self.params['N_plot'], 'Number of states to print into Spin_distribution.dat')
        input_string += self.input_line(self.params['prediag_hamiltonian'], 'Write prediagonalized Hamiltonian to PD_HAMIL.dat')
        input_string += self.input_line(self.params['eigenvectors'], 'Write eigenvectors to EIGENVEC.dat')

        return input_string

    @staticmethod
    def load_input(input_file: str):
        """Load the Hamiltonian from a file.

        Args
        -----
        input_string: str
            Hamiltonian input file
        """
        print('WARNING: Loading Hamiltonian from file is not tested')        
        
        infile = open(input_file, 'r')

        params = {}
        _ = infile.readline()
        _ = infile.readline()
        _ = infile.readline()
        
        params['Nm'] = int(infile.readline().split()[0])
        params['Spins'] = []
        
        for _ in range(params['Nm']):
            Spin = {}
            _ = infile.readline()
            Spin['S'] = float(infile.readline( ).split()[0])
            Spin['Stephen'] = [float(i) for i in infile.readline().split()[:4]]
            Spin['Stephen_axis'] = [float(i) for i in infile.readline().split()[:3]]
            Spin['H'] = [float(i) for i in infile.readline().split()[:3]]
            Spin['G'] = [float(i) for i in infile.readline().split()[:3]]
            params['Spins'].append(Spin)
        
        _ = infile.readline()
        params['Npairs'] = int(infile.readline().split()[0])
        params['pairs'] = []
        
        for _ in range(params['Npairs']):
            pair = {}
            pair['pair'] = [int(i) for i in infile.readline().split()]
            pair['J_exc'] = [float(i) for i in infile.readline().split()]
            params['pairs'].append(pair)

        _ = infile.readline()
        params['eps_QD'] = float(infile.readline().split()[0])
        params['Hubbard'] = float(infile.readline().split()[0])

        _ = infile.readline()
        params['output_file'] = infile.readline().split()[0]
        params['N_plot'] = int(infile.readline().split()[0])
        params['prediag_hamiltonian'] = F90Input.string2bool(infile.readline().split()[0])
        params['eigenvectors'] = F90Input.string2bool(infile.readline().split()[0])
        
        infile.close()
        
        return params

    def create_output_dict(self):
        
        output_dict = {}
        output_dict['hamiltonian'] = self.params['output_file']
        output_dict['eigvalues'] = 'Eigenvalues.dat'

        if self.params['eigenvectors']:
            output_dict['eigenvectors'] = 'Eigenvectors.dat'
        if self.params['prediag_hamiltonian']:
            output_dict['pre-diag_hamiltonian'] = 'PD_HAMIL.dat'
        
        return output_dict
    
    def load_output(self, output_dict):
        result_dict = {}
        result_dict['eigenvalues'] = np.loadtxt('Eigenvalues.dat')[:,1]
        return result_dict
        

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


class Floquet(F90Input): 

    def __init__(self, floquet_dict: Union[dict, str],
                 line_lenght = 80, padding_lenght = 40):
        super(Floquet, self).__init__(line_lenght, padding_lenght)

        assert isinstance(floquet_dict, (dict, str)), "Floquet input should be a dictionary or a file name"
        if isinstance(floquet_dict, str):
            assert os.path.exists(floquet_dict), f"Floquet input file {floquet_dict} does not exist"
            self.params = self.load_input(floquet_dict)
        else:
            self.params = floquet_dict

        self.check_dictionary(self.params, floq_keys, 'Floquet input')
        self.code_version = 'floquet'

        print('WARNING: No error capture for Floquet input. Worst case scenario is that the code will crash.')

    def write_input(self,):

        input_string = self.create_header('', '*')
        input_string += self.create_header('Floquet Input', '*')
        input_string += self.create_header('', '*')
        input_string += self.input_line(self.params['n_max'], 'Max order of the Fourier series')
        input_string += self.input_line(self.params['frequency'], 'Frequency of the driving field (GHz)')

        input_string += self.create_header('Electrode set-up ', '-')
        input_string += self.input_line(self.params['gamma0'][0], 'Gamma_R_0 (meV)')
        input_string += self.input_line(self.params['gamma0'][1], 'Gamma_L_0 (meV)')
        input_string += self.input_line(self.params['A'][0], 'A_R')
        input_string += self.input_line(self.params['A'][1], 'A_L')
        input_string += self.input_line(self.params['phi'], 'Phase of the driving field (radians)')
        input_string += self.input_line(self.params['seha'], 'Multiplies the 2nd Harmonic (UNUSED)')
        input_string += self.input_line(self.params['cutoff'], 'Cutoff for integral Lambshift (meV)')
        input_string += self.input_line(self.params['gammaC'], 'Broadening of Green\'s function (meV)')
        input_string += self.input_line(self.params['integral_points'], 'Number of points for I11 and I21')
        input_string += self.input_line(self.params['gwidth'], 'Width of the Gaussian for PDOS (meV)')
        input_string += self.input_line(self.params['gau'], 'Width of the Gaussian for PDOS (meV)')
        input_string += self.input_line(self.params['norb'], 'Number of orbital, 1 will not use any Diego Lehman coeff')
        input_string += self.input_line(self.params['feedback'], 'move the tip to set current Iset')
        input_string += self.input_line(self.params['Iset'], 'Set current (pA) if feedback is true')
        input_string += self.input_line(self.params['Itol'], 'Tolerance for the Iset current (pA)')
        input_string += self.input_line(self.params['ratio'], 'Ratio between gammas for the feedback')

        input_string += self.create_header('Bias, temperature, and spin polarization', '-')
        input_string += self.input_line(self.params['fermiP'], '0 EF closer to muR if gR>>gL, 1 EF to muL, 2 Ef at mu with gamma>> (UNUSED)')
        input_string += self.input_line(self.params['bias'][0], 'Right elctrode bias (mV)')
        input_string += self.input_line(self.params['bias'][1], 'Left electrode bias (mV)')
        input_string += self.input_line(self.params['Temperature'], 'Temperature (K)')
        input_string += self.input_line(self.params['Spin_polarization'][0], 'Right electrode spin polarization')
        input_string += self.input_line(self.params['Spin_polarization'][1], 'Left electrode spin polarization')
        input_string += self.input_line(self.params['Electrode'], 'Current measurement: 0 is left and 1 is right electrode')

        input_string += self.create_header('Bessel function', '-')
        input_string += self.input_line(self.params['bessel_amplitude'][0], 'B_R strengt of the time depenndet pulse for right electrode')
        input_string += self.input_line(self.params['bessel_amplitude'][1], 'B_L strengt of the time depenndet pulse for left electrode')
        input_string += self.input_line(self.params['p_max'], 'Max order of Bessel function in both directions')

        input_string += self.create_header('Output', '-')
        input_string += self.input_line(self.params['write_populations'], 'Write populations')
        input_string += self.input_line(self.params['write_coherence'], 'Write coherences')
        input_string += self.input_line(self.params['spinflo'], 'Write spin')

        input_string += self.create_header('Misc', '-')
        input_string += self.input_line(self.params['redimension'], 'Redimension the Hamiltonian')
        input_string += self.input_line(self.params['Nd'], 'Redimension number - must include empty states for non-zero current')

        return input_string
    
    @staticmethod
    def load_input(input_file: str):
        
        infile = open(input_file, 'r')
        params = {}

        _ = infile.readline()
        _ = infile.readline()
        _ = infile.readline()

        params['n_max'] = int(infile.readline().split()[0])
        params['frequency'] = float(infile.readline().split()[0])
        
        _ = infile.readline()
        params['gamma0'] = [float(infile.readline().split()[0]), 
                            float(infile.readline().split()[0])]
        
        params['A'] = [Floquet.load_complex(infile.readline()),
                        Floquet.load_complex(infile.readline())]
        params['phi'] = float(infile.readline().split()[0])
        params['seha'] = float(infile.readline().split()[0])
        params['cutoff'] = float(infile.readline().split()[0])
        params['gammaC'] = Floquet.load_complex(infile.readline())
        params['integral_points'] = int(infile.readline().split()[0])
        params['gwidth'] = float(infile.readline().split()[0])
        params['gau'] = float(infile.readline().split()[0])
        params['norb'] = int(infile.readline().split()[0])
        params['feedback'] = F90Input.string2bool(infile.readline().split()[0])
        params['Iset'] = float(infile.readline().split()[0])
        params['Itol'] = float(infile.readline().split()[0])
        params['ratio'] = float(infile.readline().split()[0])
        
        _ = infile.readline()
        params['fermiP'] = int(infile.readline().split()[0])
        params['bias'] = [float(infile.readline().split()[0]), 
                          float(infile.readline().split()[0])]
        params['Temperature'] = float(infile.readline().split()[0])
        params['Spin_polarization'] = [float(infile.readline().split()[0]), 
                                       float(infile.readline().split()[0])]
        params['Electrode'] = int(infile.readline().split()[0])

        _ = infile.readline()
        params['bessel_amplitude'] = [float(infile.readline().split()[0]),
                                        float(infile.readline().split()[0])]
        params['p_max'] = int(infile.readline().split()[0])

        _ = infile.readline()
        params['write_populations'] = F90Input.string2bool(infile.readline().split()[0])
        params['write_coherence'] = F90Input.string2bool(infile.readline().split()[0])
        params['spinflo'] = F90Input.string2bool(infile.readline().split()[0])
        
        _ = infile.readline()
        params['redimension'] = F90Input.string2bool(infile.readline().split()[0])
        params['Nd'] = int(infile.readline().split()[0])
        infile.close()

        return params

    def create_output_dict(self):
        output_dict = {}

        output_dict['hamiltonian'] = 'Hamiltonian.dat'
        output_dict['Eval'] = 'Eigenvalues.dat'
        output_dict['Evec'] = 'Eigenstates.dat'
        output_dict['Spin_dist'] = 'Spin_distribution.dat'

        output_dict['current'] = 'Current_0.dat'
        output_dict['rates'] = 'rates_floquet.dat'
        output_dict['rates_avg'] = 'rates_floquet0.dat'

        if self.params['write_populations']:
            output_dict['populations'] = 'POPULATIONS.dat'
        if self.params['write_coherence']:
            output_dict['coherence'] = 'COHERENCES.dat'
        if self.params['spinflo']:
            output_dict['spin'] = 'SpinFloquet.dat'
        
        return output_dict

    def load_output(self, output_dict):
        results_dict = {}

        I = np.loadtxt(output_dict['current'], skiprows=1)
        results_dict['DC'] = I[np.where(I[:,0] == 0),1][0,0]

        return results_dict
