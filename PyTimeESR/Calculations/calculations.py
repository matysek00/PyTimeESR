import numpy as np 
from ase.units import Hartree, _hplanck, J

from .fortran import h_qd, qme_f
from ..misc import find_bessel_order

BohrMagneton = 5.7883818066E-5 # eV/T
hplanck = _hplanck*J # eV*s
mili = 1E-3
Giga = 1E9 

def Solve_H(eps: float, U: float, 
            magnetic_field: np.array, gyro: np.array, 
            Spins: np.array = np.array([.5]), 
            spinaxis: np.array = np.array([[0,0,1]]),
            Stephen: np.array = np.array([[0.,0.,0.,0.]]),
            Jexch: np.array = np.array([[0.,0.,0.]]),
            pairs: np.array = np.array([[],[]]),
            verbose: bool = False
            ) ->  (np.array, _np.array, np.array, np.array, np.array): 
    """Runs the Hamiltonia of the Unpertrubed QM Dot. 

    Parameters:
    ----------
    eps: float
        Single electron energy splitting (meV) 
    U: float
        Coulomb repulsion energy (meV)
    magnetic_field: np.array(Nmolecules,3) (Tesla)
        Magnetic field vector on each molecule
    gyro: np.array(Nmolecules,3)
        Gyromagnetic ratio vector for each molecule (meV)
    Spins: np.array(Nmolecules,)
        Array of spin values for each molecule, default is [0.5] (single molecule), 
        Spins[0] must be 0.5 for transport. 
    spinaxis: np.array(Nmolecules,3)
        Array of spin axis for each molecule, default is [[0,0,1]] (z-axis, single molecule)
    Stephen: np.array(Nmolecules,4)
        Stevens parameters for each molecule [B20, B22, B40, B44], default is [[0.,0.,0.,0.]] (no anisotropy)
    Jexch: np.array(Npairs,3)
        Exchange interaction matrix, default is [] (no interaction) (meV*GHz)
    pairs: np.array(2,Npairs)
        Array of index pairs for exchange interaction, default is [[],[]] (no interactions)
    verbose: bool
        If True, prints additional information. 
    
    Returns:
    -------
    lamb: np.array(N,N,2)
        Overlap matrix between spin and impurity states
        lamda(i,j,0) = <i|S+|j>, lamda(i,j,1) = <i|S-|j> (double check)
    Delta: np.array(N,N)
        Delta(i,j) = Eigen(i) - Eigen(j)
    H: np.array(N,N)
        Hamiltonian matrix
    Eigen: np.array(N)
        Eigenvalues of the Hamiltonian 
    Ss: np.array(Nmol, 3, N, N)
        Spin matrix elements for each molecule
    """

    # check dimensions 
    Nmolecules = len(Spins)
    Npairs = pairs.shape[1]
    
    if verbose:
        print("Solving Hamiltonian for", Nmolecules, "molecules with", Npairs, "exchange pairs.")
    
    assert magnetic_field.shape == (Nmolecules,3), f'magnetic_field must be of shape (Nmolecules={Nmolecules:d},3)'
    assert gyro.shape == (Nmolecules,3), f'gyro must be of shape (Nmolecules={Nmolecules:d},3)'
    assert spinaxis.shape == (Nmolecules,3), f'spinaxis must be of shape (Nmolecules={Nmolecules:d},3)'
    assert Stephen.shape == (Nmolecules,4), f'Stephen must be of shape (Nmolecules={Nmolecules:d},4)'
    assert Jexch.shape == (Npairs,3), f'Jexch must be of shape (Npairs={Npairs:d},3)'
    assert Spins[0] == 0.5, 'First spin must be 0.5 for transport calculations.'

    #TODO: move conversions to be external
    eps *= mili/Hartree # convert meV to Hartree
    U *= mili/Hartree # convert meV to Hartree
    
    magnetic_field *= gyro*BohrMagneton # eV 
    magnetic_field /= Hartree # convert eV to Hartree

    Stephen /= Hartree # convert meV to Hartree

    Jexch /= hplanck*Giga/Hartree # convert GHz to Hartree

    N, H_el, Nin,  Nblock, = h_qd.create_basis(eps, U, Spins)

    if verbose: 
        print('')
        print('Total basis size N:', N)
        print('Key energies for transport:')
        print('Ionization energy:',eps*1000*Hartree, 'meV')
        print('Coulomb repulsion:', U*1000*Hartree, 'meV')  
        print('On site energy plus Ionization:', (U+eps)*1000*Hartree, 'meV')
        print('')
    

    B = Stephen.T
    
    lamb, Delta, H, Eigen, Ss = h_qd.hamiltonian(
        N, Nin, Nblock, H_el, magnetic_field,spinaxis, Spins,
        B[0], B[1], B[2], B[3], Jexch, pairs[0], pairs[1])
    
    Delta *= Hartree/mili # convert back to meV
    Eigen *= Hartree/mili # convert back to meV
    H *= Hartree/mili # convert back to meV

    return lamb, Delta, H, Eigen, Ss


def Calculate_rates(frequency: float, gamma0: float,  lamb: np.array,  
    spin_pol: float, bias: float,Delta: np.array, Adrive: float=0, 
    Vrf: float=0., transport_exponent: float = 0.0, Temperature: float = .5, 
    Nfour: int = 3, broadening: float = 1E-3, cutoff: float = 500., Nint: int = 1000,
    verbose: bool = False,)-> (np.array, np.array):
    """Calculates the transition rates between states due to tunneling and driving fields.

    Parameters:
    ----------
    frequency: float
        Driving frequency (GHz)
    gamma0: float
        Base tunneling rate (meV)
    lamb: np.array(N,N,2)
        Overlap matrix between spin and impurity states
        lamda(i,j,0) = <i|S+|j>, lamda(i,j,1) = <i|S-|j>
    spin_pol: float
        Spin polarization of the leads (0 to 1)
    bias: float
        Bias voltage (meV)
    Delta: np.array(N,N)
        Delta(i,j) = Eigen(i) - Eigen(j)
    Adrive: float
        Amplitude of the driving field, default is 0 
    Vrf: float
        RF voltage amplitude, default is 0 (meV)
    transport_exponent: float
        Exponent for transport rate dependence on bias, default is 0.0 (meV^-1)
    Temperature: float
        Temperature, default is 0.5 (K)
    Nfour: int
        Number of Fourier components to consider, default is 3
    broadening: float
        Broadening of energy levels, default is 1E-3 (meV)
    cutoff: float
        Cutoff energy for integration, default is 500. (meV)
    Nint: int
        Number of integration points, default is 1000
    verbose: bool
        If True, prints additional information.
    
    Returns:
    -------
    G: np.array(N,N,N,N,2)
        Transition rates 
    G_bar: np.array(N,N,N,N,2)
        Complementary transition rates 
    """

    gamma0 *= mili/Hartree # convert meV to Hartree
    broadening *= mili/Hartree # convert meV to Hartree
    cutoff *= mili/Hartree # convert meV to Hartree
    bias *= mili/Hartree # convert meV to Hartree
    Vrf *= mili/Hartree # convert meV to Hartree
    Delta *= mili/Hartree # convert meV to Hartree
    frequency *= hplanck/Hartree*Gig # convert GHz to Hz

    
    pmax = find_bessel_order(Vrf, frequency)
    
    if verbose:
        print(f'Using pmax={pmax:d} for Bessel function expansion with Vrf={Vrf:.3f} meV and frequency={frequency:.3f} GHz.')
    
    NF = 2*Nfour + 1
    G, Gbar = qme_f.rates(frequency, gamma0, lamb, spin_pol, NF, pmax, 
        Adrive, Vrf, broadening, bias, Delta, cutoff, Temperature,
        Nint, transport_exponent)


    return G, Gbar