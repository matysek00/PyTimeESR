# PyTimeESR

A simple pythonic wrapper around the TimeESR code

## Minimal example

```
code_path = 'path_to_TimeESR_src'
run_path = 'path_to_run'

fn_input_qd = 'path_to_QD.input'
fn_input_dyn = 'path_to_TimeESR.input'

# make 
PyTimeESR.make(code_path)

# initialize Simulation
Sim = PyTimeESR.Simulation(fn_input_qd, fn_input_dyn, run_path, code_path)

# change parameters from the file
# see PyTimeESR.default_ham and PytimeESR.default_dyn for params structure 
Sim.Dynamics.params['intervals'][0]['freq'][0]['Frequency'] = 17.0 #GHz
Sim.Dynamics.params['intervals'][0]['freq'][0]['Amplitude'] = .5 

# run the simulation
Sim.run()
```

To use the Floquet version of the code, we use `version`
```
fn_input_floq = 'path_to_floquet.input'
Sim = PyTimeESR.Simulation(fn_input_ham, fn_input_floq, run_path, 
            code_path, code_version='floquet')
Sim.run()
```
Alternatively one can use dictionaries instead of input files, which might be preferable if you want edit multiple parameters. 
```
ham_dict = PyTimeESR.Hamiltonian.load_input(fn_ham)
floq_dict = PyTimeESR.Floquet.load_input(fn_dyn)

floq_dict['n_max'] = 150
floq_dict['p_max'] = 210
floq_dict['bessel_amplitude'] = [0., 6.5]
Sim = PyTimeESR.Simulation(ham_dict, floq_dict, run_path, 
            code_path, code_version='floquet')
```
For Bessel version of the TimeESR code use `version='bessel'`

```
ham_dict = PyTimeESR.Hamiltonian.load_input(fn_ham)
dyn_dict = PyTimeESR.Floquet.load_input(fn_dyn, version='bessel')

Sim = PyTimeESR.Simulation(ham_dict, dyn_dict, run_path, 
            code_path, code_version='bessel')
```

There are some potentially useful functions in `PyTimeESR.empirical` but this pretty much everything. 
I put some more examples, everything should be in this file already.