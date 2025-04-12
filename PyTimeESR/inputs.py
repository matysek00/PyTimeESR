
from .default_inputs import *
               
class F90Input(): 
    """Parent class to write fortran input.
    """
    
    fmt = {
        tfloat: '{:.6f}',
        tint: '{:d}',
        tstr: '{}',
        tbool: '{}',
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
        ttype = tfloat if t in tfloat else tint if t in tint else tstr if t in tstr else None
        assert ttype is not None, f"Type {t} not supported in input line"
        
        line = ' '.join([self.fmt[ttype].format(i) for i in value])
        pad = self.padding_lenght - len(line)
        pad = 0 if pad < 0 else pad
        line += ' ' * pad + '! ' + comment + '\n'
        return line
    

class Hamiltonian(F90Input):
    """Write input for the Hamiltonian. 

    Args
    -----
    hamiltonian_dict: dict
        Dictionary with the Hamiltonian input. 
    """
    
    def __init__(self, hamiltonian_dict: dict, line_lenght = 80, padding_lenght = 40):
        super(Hamiltonian, self).__init__(line_lenght, padding_lenght)
        self.check_dictionary(hamiltonian_dict, ham_keys, 'Hamiltonian input')
        
        Nm = hamiltonian_dict['Nm']
        Npairs = hamiltonian_dict['Npairs']
        
        assert  Nm > 0, "Number of spins should be greater than 0"
        assert Npairs >= 0, "Number of pairs should be greater than or equal to 0"
        
        assert len(hamiltonian_dict['Spins']) == Nm, "Number of spins should match the number of spins in the dictionary"
        assert len(hamiltonian_dict['pairs']) == Npairs, "Number of pairs should match the number of pairs in the dictionary"

        for i in range(Nm): 
            self.check_dictionary(hamiltonian_dict['Spins'][i], ham_keys_spins, 'Hamiltonian input - Spins')
        for i in range(Npairs): 
            self.check_dictionary(hamiltonian_dict['pairs'][i], ham_keys_pairs, 'Hamiltonian input - Pairs')
        
        S_trans = hamiltonian_dict['Spins'][0]['S']
        assert S_trans == .5, 'Transport electron should be a spin 1/2'
        
        self.params = hamiltonian_dict

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


class Dynamics(F90Input):
    """Write input for the dynamics. 

    Args
    -----
    dynamics_dict: dict
        Dictionary with the dynamics input. 
        """

    def __init__(self, dynamics_dict: dict, line_lenght = 80, padding_lenght = 40):
        super(Dynamics, self).__init__(line_lenght, padding_lenght)
        self.check_dictionary(dynamics_dict, dyn_keys, 'Dynamics input')

        N_interval = dynamics_dict['N_interval']
        Nfreq = dynamics_dict['Nfreq']
        Nbias = dynamics_dict['Nbias']
    
        assert N_interval > 0, "Number of intervals should be greater than 0"
        assert Nfreq > 0, "Number of frequencies should be greater than 0"
        
        assert len(dynamics_dict['intervals']) == N_interval, "Number of intervals should match the number of intervals in the dictionary"
        assert len(dynamics_dict['biases']) == Nbias, "Number of biases should match the number of biases in the dictionary"

        for i in range(N_interval): 
            self.check_dictionary(dynamics_dict['intervals'][i], dyn_keys_interval, 'Dynamics input - Intervals')
            assert len(dynamics_dict['intervals'][i]['freq']) == Nfreq, "Number of frequencies should match the number of frequencies in the dictionary"
            for j in range(Nfreq):
                self.check_dictionary(dynamics_dict['intervals'][i]['freq'][j], dyn_keys_freq, 'Dynamics input - Frequencies')
        for i in range(Nbias):
            self.check_dictionary(dynamics_dict['biases'][i], dyn_keys_bias, 'Dynamics input - Biases')

        intervals = dynamics_dict['intervals']
        assert all(intervals[i]['t0'] < intervals[i]['tf'] for i in range(N_interval)), "Interval times should be in increasing order"
        assert all(intervals[i]['t0'] == intervals[i-1]['tf'] for i in range(1, N_interval)), "Interval times should be continuous"
        assert intervals[0]['t0'] == 0, "First interval should start at 0"
        assert intervals[-1]['tf'] == dynamics_dict['t_final'], "Last interval should end at t_final"
        
        if  dynamics_dict['use_bessel']:
            assert Nfreq == 1, "Number of frequencies should be 1 for Bessel function"
            assert N_interval == 1, "Number of intervals should be 1 for Bessel function"
            assert dynamics_dict['p_max'] > dynamics_dict['n_max'], "Max order of Bessel function should be greater than max frequency"
            assert dynamics_dict['n_max'] >= 0, "Max frequency of Bessel function should be greater than or equal to 0"

        self.params = dynamics_dict

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

        input_string += self.create_header('Bessel function', '-')
        input_string += self.input_line(self.params['use_bessel'], 'Use Bessel function')
        input_string += self.input_line(self.params['bessel_aplitude'][0], 'B_R strengt of the time depenndet pulse for right electrode')
        input_string += self.input_line(self.params['bessel_aplitude'][1], 'B_L strengt of the time depenndet pulse for left electrode')
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
