
import os
import sys

from . import inputs
class Simulation(): 
    """Create, run, and analyze TimeESR simulation.

    Args
    -----
    """
    output_dict = {
        'population': 'POPULATIONS.dat',
        'spin_dynamics': 'SpinDynamics.dat',
        'esr': 'ESR.dat',
        'hamiltonian': 'Hamiltonian.output',
        'pre-diag_hamiltonian': 'PD_HAMIL.dat',
        'spin_distribution': 'Spin_distrubution.dat',
        'current': 'Current.dat',
        'population_average': 'POP_AVE.dat',
    }

    def __init__(self, Ham_dict: dict, Dyn_dict: dict, run_path: str, code_path: str):
        
        """Initialize the simulation with Hamiltonian and dynamics dictionaries.

        Args
        -----
        Ham_dict : dict
            Dictionary containing Hamiltonian parameters.
        Dyn_dict : dict
            Dictionary containing dynamics parameters.
        """
        
        self.Ham = inputs.Hamiltonian(Ham_dict)
        self.Dyn = inputs.Dynamics(Dyn_dict)


        exec_path = os.path.join(code_path, 'TimeESR.x')

        assert os.path.exists(run_path), f"Run path {run_path} does not exist."
        assert os.path.exists(exec_path), f"Executable path {exec_path} does not exist."
        
        self.run_path = run_path
        self.exec_path = exec_path

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
        os.system(self.exec_path)
        os.chdir(current_path)


def make(path): 
    """Compile the TimeESR code.
    """
    pwd = os.getcwd()
    os.chdir(path)
    os.system('make clean')
    os.system('make')
    os.chdir(pwd)