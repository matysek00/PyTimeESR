import numpy as np

# Define some general types
tfloat = (float, np.float32, np.float64)
tcomplex = (complex, np.complex64, np.complex128)
tint = (int, np.int32, np.int64)
tlist = (list, np.ndarray)
tstr = (str, np.str_)
tbool = (bool, np.bool_)

# Keys and thier types, if the type is a list, the first element
# is the type of the list, the second element type of elememts, 
# and the third element is the number of elements in the list.
ham_keys = {
    'Nm': tint, 'Spins': [list, dict, None], 'Npairs': tint, 'pairs': [list, dict, None],
    'eps_QD': tfloat, 'Hubbard': tfloat,
    'output_file': str, 'N_plot': tint, 'prediag_hamiltonian': tbool,
    'eigenvectors': tbool}
ham_keys_spins = {
    'S': tfloat, 'Stephen': [list, tfloat, 4], 'Stephen_axis': [list, tfloat, 3],
    'H': [list, tfloat, 3], 'G': [list, tfloat, 3]}   
ham_keys_pairs = {
    'pair': [list, tint, 2], 'J_exc': [list, tfloat, 3],}

default_ham = {
### Number of spins
    'Nm': 1,
### Spin Properties
    'Spins': [{
        'S': .5,
        'Stephen': [0., 0., 0., 0.],
        'Stephen_axis': [0., 0., 1.],
        'H': [.60827625,0.,0.],
        'G': [2., 2., .2]}],
### Exchnge 
    'Npairs': 0, 
    'pairs': [
#        {'pair': [0,1],
#        'J_exc': [0.0, 0.0, 0.0],},
        ],
### Electronic interactions
    'eps_QD': -10.,
    'Hubbard': 100., 
### Output 
    'output_file': 'Hamiltonian.output',
    'N_plot': 4,
    'prediag_hamiltonian': False,
    'eigenvectors': False,
}

dyn_keys = {
    'Ntime': tint, 't_initial': tfloat, 't_final': tfloat,
    'N_interval': tint, 'Nfreq': tint, 'intervals': [list, dict, None],
    'gamma0': [list, tfloat, 2], 'gamma1': [list, tfloat, 2],
    'cutoff': tfloat, 'gammaC': tfloat, 'integral_points': tint,
    'Nbias': tint, 'biases': [list, dict, None], 'Temperature': tfloat,
    'spin_polarization': [list, tfloat, 2], 'Electrode': tint,
    'use_bessel': tbool, 'bessel_amplitude': [list, tfloat, 2],
    'p_max': tint, 'n_max': tint,
    'population': tbool, 'density_matrix': tbool,
    'output_file': str, 'output_fourier': str,
    'output_ESR': str, 'runs': tbool,
    'spindyn': tbool, 'redimension': tbool,
    'Nd': tint,
}

dyn_keys_interval = {
    't0': tfloat, 'tf': tfloat, 'freq': [list, dict, None], 'Phase': tfloat,
}
dyn_keys_freq = {
    'Amplitude': tfloat, 'Frequency': tfloat, 
}   
dyn_keys_bias = {
    'bias': [list, tfloat, 2], 'b_time': tfloat,
}

default_dyn = {
### Total Time
    'Ntime': 500000,
    't_initial': 0.0,
    't_final': 150.0,
### Pulses
    'N_interval': 1, 
    'Nfreq': 1,
    'intervals': [{
        't0': 0.0,
        'tf': 150.0,
        'freq': [{
            'Amplitude': 1.0,
            'Frequency': 16.0,
        },],
        'Phase': 0.0,
        },],
### Electrodes [R,L]
    'gamma0': [0.0020, 0.0001],
    'gamma1': [0.0000, 0.0002],
    'cutoff': 100., 
    'gammaC': 0.0020,
    'integral_points': 500000,
### bias, temperature, spin polarization
    'Nbias': 1,
    'biases': [{
        'bias': [15., -15.],
        'b_time': 150.
        },],
    'Temperature': .5,
    'spin_polarization': [0., 1.],
    'Electrode' : 0, 
### Bessel functions
    'use_bessel': False,
    'bessel_amplitude': [0.0, 0.0],
    'p_max': 10,
    'n_max': 3,
### Output
    'population': False,
    'density_matrix': False, 
    'output_file': 'C.dat',
    'output_fourier': 'S.dat',
    'output_ESR': 'ESR.dat',
### Read previous hamiltonian
    'runs': False, 
### Misc
    'spindyn': False, 
    'redimension': False,
    'Nd': 4,
    }

floq_keys = {
    'n_max': tint, 'frequency': tfloat,
    'gamma0': [list, tfloat, 2], 'A': [list, tcomplex, 2],
    'phi': tfloat, 'seha': tfloat, 'cutoff': tfloat, 'gammaC': tcomplex,
    'integral_points': tint, 'gwidth': tfloat, 'gau': tfloat, 'norb': tint,
    'feedback': tbool, 'Iset': tfloat, 'Itol': tfloat, 'ratio' : tfloat,
    'fermiP' : tint, 'bias': [list, tfloat, 2], 'Temperature': tfloat,
    'Spin_polarization': [list, tfloat, 2], 'Electrode': tint,
    'bessel_amplitude' : [list, tfloat, 2], 'p_max': tint, 
    'write_populations': tbool, 'write_coherence': tbool,
    'spinflo': tbool, 'redimension': tbool, 'Nd': tint,
}

default_floq = {
    'n_max': 3,
    'frequency': 16.0,
    'gamma0': [0.0020, 0.0001],
    'A': [complex(0.0), complex(0.0002)],
    'phi': 0.0,
    'seha': 1.,
    'cutoff': 100.,
    'gammaC': complex(0.0020),
    'integral_points': 500000,
    'gwidth': 100.,
    'gau': 0.0,
    'norb': 1,
    'feedback': False,
    'Iset': 1.0,
    'Itol': 0.001,
    'ratio': 10.,
    'fermiP': 0,
    'bias': [15., -15.],
    'Temperature': 0.5,
    'Spin_polarization': [0., 1.],
    'Electrode': 0,
    'bessel_amplitude': [0.0, 0.0],
    'p_max': 10,
    'write_populations': False,
    'write_coherence': False,
    'spinflo': False,
    'redimension': False,
    'Nd': 4,
}