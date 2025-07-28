
import os
import sys
import time

import numpy as np

from . import inputs

class Simulation(): 
    """Create, run, and analyze TimeESR simulation.

    Args
    -----
    """
    output_dict = {
        'spin_distribution': 'Spin_distrubution.dat',
        'current': 'Current.dat',
        'population_average': 'POP_AVE.dat'}

    results_dict = {
        'esr': None, 
        'run_time': None,}

    def __init__(self, Ham_dict: dict, Dyn_dict: dict, 
                 run_path: str, code_path: str, 
                 code_version: str = 'standart'):
        
        """Initialize the simulation with Hamiltonian and dynamics dictionaries.

        Args
        -----
        Ham_dict : dict
            Dictionary containing Hamiltonian parameters.
        Dyn_dict : dict
            Dictionary containing dynamics parameters. Or the Floquet parameters if code_version is 'floquet'.
        run_path : str
            Path to the directory where the simulation will be run.
        code_path : str
            Path to the TimeESR code directory.
        code_version : str
            Version of the TimeESR code to use 'standart' (default) or 'bessel'.
        """
        
        assert code_version in ['bessel', 'standart', 'floquet'], \
            f"Code version {code_version} is not supported. Use 'bessel' or 'standart'."

        executable = 'Floquet_ESR_v7.4.0.out' if code_version == 'floquet' else 'TimeESR.x'
        exec_path = os.path.join(code_path, executable)

        assert os.path.exists(run_path), f"Run path {run_path} does not exist."
        assert os.path.exists(exec_path), f"Executable path {exec_path} does not exist."
        
        self.run_path = run_path
        self.exec_path = exec_path

        self.Ham = inputs.Hamiltonian(Ham_dict)
        if code_version == 'floquet':
            self.Dyn = inputs.Floquet(Dyn_dict)
        else: 
            self.Dyn = inputs.Dynamics(Dyn_dict, code_version=code_version)

        self.output_dict = {**self.output_dict,
                            **self.Dyn.create_output_dict(), 
                            **self.Ham.create_output_dict()} 
        
        self.output_dict = {key: os.path.join(run_path, value) for key, 
                            value in self.output_dict.items()}

    def run(self, outfile = None):
        """Run the TimeESR simulation.
        """

        current_path = os.getcwd()
        
        dyn_fn = 'Floquet.input' if self.Dyn.code_version == 'floquet' else 'TimeESR.input'
        fnham = os.path.join(self.run_path, 'H_QD.input')
        fnesr = os.path.join(self.run_path, dyn_fn)

        fham = open(fnham, 'w')
        fham.write(self.Ham.write_input())
        fham.close()

        fesr = open(fnesr, 'w')
        fesr.write(self.Dyn.write_input())
        fesr.close()

        command = f'{self.exec_path}'
        if outfile is not None: 
            command += f' >> {outfile}'
        
        os.chdir(self.run_path)
        t1 = time.time()
        os.system(f'{self.exec_path} >> {outfile}')
        t2 = time.time()
        os.chdir(current_path)

        self.results_dict['run_time'] = t2 - t1
        #self.load_output()

    def load_output(self):
        self.results_dict = {**self.results_dict,
                             **self.Ham.load_output(self.output_dict), 
                             **self.Dyn.load_output()}
    
    def get_fidelity(self, phi):
        """Calculate the fidelity between the current and a reference state.

        Args
        -----
        phi : np.ndarray
            The reference state.

        Returns
        -------
        float
            The fidelity between the current and reference states.
        """
        assert 'population' in self.output_dict, \
            'Population data not found in output dictionary.'
        
        fnpop = self.output_dict['population']
        time, F = self.fidelity_evolution(phi, fnpop)
        
        self.results_dict['fidelity'] = F
        self.results_dict['time'] = time

    def get_entropy(self,):
        """Calculate the entropy of the system.

        Returns
        -------
        float
            The entropy of the system.
        """
        assert 'population' in self.output_dict, \
            'Population data not found in output dictionary.'
        
        fnpop = self.output_dict['population']
        time, S = self.entropy_evolution(fnpop)
        
        self.results_dict['entropy'] = S
        self.results_dict['time'] = time


def make(path): 
    """Compile the TimeESR code.
    """
    pwd = os.getcwd()
    os.chdir(path)
    os.system('make clean')
    os.system('make')
    os.chdir(pwd)