import os

from typing import Union

from .F90_input import F90Input
from .default_inputs import *
from ..misc import find_bessel_order


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
        self.output_dict = self.create_output_dict()


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
        input_string += self.input_line(self.params['cutoff'], 'Cutoff for integral Lambshift (meV)')
        input_string += self.input_line(self.params['gammaC'], 'Broadening of Green\'s function (meV)')
        input_string += self.input_line(self.params['integral_points'], 'Number of points for I11 and I21')
        input_string += self.input_line(self.params['norb'], 'Number of orbital, 1 will not use any Diego Lehman coeff')
        
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
        input_string += self.input_line(self.params['p_max'][0], 'Max order of Bessel function in both directions, right')
        input_string += self.input_line(self.params['p_max'][1], 'Max order of Bessel function in both directions, left')

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
        params['cutoff'] = float(infile.readline().split()[0])
        params['gammaC'] = Floquet.load_complex(infile.readline())
        params['integral_points'] = int(infile.readline().split()[0])
        params['norb'] = int(infile.readline().split()[0])
        
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
        params['p_max'] = [int(infile.readline().split()[0]),
                           int(infile.readline().split()[0])]

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

    def load_output(self):
        results_dict = {}

        I = np.loadtxt(self.output_dict['current'], skiprows=1)
        
        if I.size == 2: # In case that only DC current was taken
            results_dict['DC'] = I[1]
        else:
            results_dict['DC'] = I[0,1]

        return results_dict
    
    def load_rates(self, output_dict, ndim):
        fn_rates = output_dict['rates']
        nmax = self.params['n_max'] 
        nfour = 2*nmax + 1 

        rates_discrete = np.loadtxt(fn_rates, skiprows=1)
        ratesR = np.zeros((ndim, ndim, ndim, ndim, nfour), dtype = np.complex128)
        ratesL = ratesR.copy()             

        for r in rates_discrete:
            ind = r[:5].astype(int)
            ind[4] += nmax
            ind[:4] -= 1

            ratesL[tuple(ind)] = np.complex128(*r[5:7])
            ratesR[tuple(ind)] = np.complex128(*r[7:9])
        
        return ratesL, ratesR
    
    def set_pmax(self):

        Vrf = self.params['bessel_amplitude'] 
        omega = self.params['frequency']

        self.params['p_max'][0] = find_bessel_order(Vrf[0], omega)
        self.params['p_max'][1] = find_bessel_order(Vrf[1], omega)
    
    def Iset(self, targetI):
        pass