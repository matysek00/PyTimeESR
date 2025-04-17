
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
                 code_version: str = 'bessel'):
        
        """Initialize the simulation with Hamiltonian and dynamics dictionaries.

        Args
        -----
        Ham_dict : dict
            Dictionary containing Hamiltonian parameters.
        Dyn_dict : dict
            Dictionary containing dynamics parameters.
        run_path : str
            Path to the directory where the simulation will be run.
        code_path : str
            Path to the TimeESR code directory.
        code_version : str
            Version of the TimeESR code to use 'bessel' (default) or 'standart'.
        """
        
        exec_path = os.path.join(code_path, 'TimeESR.x')

        assert code_version in ['bessel', 'standard'], \
            f"Code version {code_version} is not supported. Use 'bessel' or 'standard'."

        assert os.path.exists(run_path), f"Run path {run_path} does not exist."
        assert os.path.exists(exec_path), f"Executable path {exec_path} does not exist."
        
        self.run_path = run_path
        self.exec_path = exec_path

        self.Ham = inputs.Hamiltonian(Ham_dict)
        self.Dyn = inputs.Dynamics(Dyn_dict, code_version=code_version)

        self.output_dict = {**self.output_dict,
                            **self.Dyn.create_output_dict(), 
                            **self.Ham.create_output_dict()} 
        
        self.output_dict = {key: os.path.join(run_path, value) for key, 
                            value in self.output_dict.items()}

    def run(self):
        """Run the TimeESR simulation.
        """

        current_path = os.getcwd()
        
        fnham = os.path.join(self.run_path, 'H_QD.input')
        fnesr = os.path.join(self.run_path, 'TimeESR.input')

        fham = open(fnham, 'w')
        fham.write(self.Ham.write_input())
        fham.close()

        fesr = open(fnesr, 'w')
        fesr.write(self.Dyn.write_input())
        fesr.close()
        
        os.chdir(self.run_path)
        t1 = time.time()
        os.system(self.exec_path)
        t2 = time.time()
        os.chdir(current_path)

        self.results_dict['run_time'] = t2 - t1
        self.load_output()

    def load_output(self):
        self.results_dict['DC'] = np.loadtxt(self.output_dict['DC'])
        

def make(path): 
    """Compile the TimeESR code.
    """
    pwd = os.getcwd()
    os.chdir(path)
    os.system('make clean')
    os.system('make')
    os.chdir(pwd)