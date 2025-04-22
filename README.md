# PyTimeESR

A simple pythonic wrapper around the TimeESR code

## Aim

* Compare different versions of the code numerically 
* Run simulations with varying parameters

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

# load outputs
DC = Sim.results_dict['DC']
print('DC: {DC}')

```