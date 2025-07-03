import os

import numpy as np
from typing import Union

from .F90_input import F90Input
from .default_inputs import *

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
